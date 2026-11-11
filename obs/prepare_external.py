"""Prepare free legs for cases whose ligands aren't in the openfe-benchmarks data,
building the ligands from each case's transformation JSON, with the partial charges
used by OpenFE. The ligands are parameterised, merged and solvated as by
openfe-benchmarks' obs/prepare_rbfe.py, once with OpenFE's mapping and once with
BioSimSpace's.

    python obs/prepare_external.py --benchmarks ../openfe-benchmarks
"""

import argparse
import io
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np


def components(path):
    """Return the two ligands, as RDKit molecules, and OpenFE's mapping."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    objects = {key: value for key, value in json.loads(Path(path).read_text())}
    mapping = next(
        v for v in objects.values() if v.get("__qualname__") == "LigandAtomMapping"
    )

    ligands = []
    for side, default in (("componentA", "A"), ("componentB", "B")):
        component = objects[mapping[side][":gufe-key:"]]
        mol = Chem.RWMol()
        charges = []
        for z, isotope, formal, aromatic, chiral, _, props, _ in component["atoms"]:
            atom = Chem.Atom(z)
            atom.SetIsotope(isotope)
            atom.SetFormalCharge(formal)
            atom.SetIsAromatic(aromatic)
            atom.SetChiralTag(Chem.ChiralType.values[chiral])
            atom.SetNoImplicit(True)
            mol.AddAtom(atom)
            charges.append(props.get("PartialCharge"))
        for i, j, order, _, _ in component["bonds"]:
            mol.AddBond(i, j, Chem.BondType.values[order])
            if order == 12:
                mol.GetBondBetweenAtoms(i, j).SetIsAromatic(True)
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)

        positions = np.load(io.BytesIO(component["conformer"][0].encode("latin-1")))
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(positions):
            conformer.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
        mol.AddConformer(conformer, assignId=True)
        Chem.AssignStereochemistryFrom3D(mol)

        name = component["molprops"].get("ofe-name", default)
        mol.SetProp("_Name", name)
        if None not in charges:
            for atom, charge in zip(mol.GetAtoms(), charges):
                atom.SetDoubleProp("PartialCharge", charge)
        ligands.append(mol)

    if ligands[0].GetProp("_Name") == ligands[1].GetProp("_Name"):
        ligands[0].SetProp("_Name", "A")
        ligands[1].SetProp("_Name", "B")
    mapping = {int(k): v for k, v in mapping["componentA_to_componentB"].items()}
    return ligands, mapping


def parameterise_with_charges(mol, BSS, _common):
    """Parameterise a molecule without stored charges, keeping the force field's."""
    from rdkit import Chem

    with tempfile.TemporaryDirectory() as tmp:
        sdf = Path(tmp) / "mol.sdf"
        with Chem.SDWriter(str(sdf)) as writer:
            writer.write(mol)
        parameterised = BSS.Parameters.parameterise(
            BSS.IO.readMolecules(str(sdf))[0], _common.DEFAULT_LIGAND_FORCEFIELD, work_dir=tmp
        ).getMolecule()
    cursor = parameterised._sire_object.cursor()
    for residue in cursor.residues():
        residue.name = _common.LIGAND_RESNAME
    parameterised._sire_object = cursor.molecule().commit()
    return parameterised


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cases", default=Path(__file__).parent / "cases_external.json")
    parser.add_argument("--benchmarks", required=True, help="openfe-benchmarks directory.")
    parser.add_argument("--output-dir", default="obs_output/external")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    sys.path.insert(0, str(Path(args.benchmarks).resolve() / "obs"))
    import BioSimSpace as BSS
    import _common
    import prepare_rbfe

    _common.setup_logging()
    root = Path(args.cases).resolve().parent.parent

    for case in json.loads(Path(args.cases).read_text()):
        path = next((root / case["case"]).glob(case["json"]))
        ligands, mapping = components(path)
        names = [mol.GetProp("_Name") for mol in ligands]

        param_dir = Path(args.output_dir) / case["id"] / "parameterised"
        (param_dir / "ligands").mkdir(parents=True, exist_ok=True)
        for mol in ligands:
            out = param_dir / "ligands" / f"{_common.safe_name(mol.GetProp('_Name'))}.bss"
            if args.overwrite or not out.exists():
                if mol.GetAtomWithIdx(0).HasProp("PartialCharge"):
                    parameterised = _common.parameterise_small_molecule(
                        mol, _common.DEFAULT_LIGAND_FORCEFIELD, _common.LIGAND_RESNAME
                    )
                else:
                    print(f"{case['id']}: no partial charges for "
                          f"{mol.GetProp('_Name')}, using the force field's", flush=True)
                    parameterised = parameterise_with_charges(mol, BSS, _common)
                BSS.Stream.save(parameterised.toSystem(), str(out.with_suffix("")))

        edge = SimpleNamespace(
            name=f"{names[0]}~{names[1]}",
            componentA=SimpleNamespace(name=names[0]),
            componentB=SimpleNamespace(name=names[1]),
            componentA_to_componentB={str(k): v for k, v in mapping.items()},
            annotations={},
        )

        for method in ("openfe", "biosimspace"):
            out = Path(args.output_dir) / case["id"] / method
            if (out / "free.bss").exists() and not args.overwrite:
                continue
            out.mkdir(parents=True, exist_ok=True)
            merged, lig_a, lig_b, used, details = prepare_rbfe.merge_edge(
                edge, method, param_dir, False, False, 1000
            )
            with tempfile.TemporaryDirectory() as tmp:
                system = _common.solvate(
                    merged.toSystem(), "1.5 nm", "dodecahedron", 0.15, work_dir=tmp
                )
            BSS.Stream.save(system, str(out / "free"))
            _common.write_json(
                out / "edge.json",
                {
                    "case": case["case"],
                    "ligand_a": names[0],
                    "ligand_b": names[1],
                    "charge_a": _common.net_charge(lig_a),
                    "charge_b": _common.net_charge(lig_b),
                    "mapping": used,
                    "mapping_method": details,
                    "openfe_charges": all(
                        m.GetAtomWithIdx(0).HasProp("PartialCharge") for m in ligands
                    ),
                },
            )
            print(f"{case['id']} {method}: {len(used)} mapped atoms", flush=True)


if __name__ == "__main__":
    main()
