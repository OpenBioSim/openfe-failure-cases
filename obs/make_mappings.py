"""Write mappings files for the cases in cases.json, for openfe-benchmarks'
obs/prepare_rbfe.py --mappings. OpenFE's atom mapping is taken from each case's
transformation JSON and checked against the benchmark ligands. Two files are
written: one that uses OpenFE's mapping, and one that maps the atoms with
BioSimSpace, falling back to OpenFE's mapping if BioSimSpace's can't be merged.

    python obs/make_mappings.py --benchmarks ../openfe-benchmarks
"""

import argparse
import json
from pathlib import Path


def openfe_mapping(case):
    """Return the ligand names and mapping from a case's transformation JSON."""
    path = next(case.glob("in_solvent_*.json"), None) or next(case.glob("in_complex_*.json"))
    objects = {key: value for key, value in json.loads(path.read_text())}
    mapping = next(v for v in objects.values() if v.get("__qualname__") == "LigandAtomMapping")
    names = []
    elements = []
    for side in ("componentA", "componentB"):
        component = objects[mapping[side][":gufe-key:"]]
        names.append(component["molprops"]["ofe-name"])
        elements.append([atom[0] for atom in component["atoms"]])
    return names, elements, {int(k): v for k, v in mapping["componentA_to_componentB"].items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cases", default=Path(__file__).parent / "cases.json")
    parser.add_argument("--benchmarks", required=True, help="openfe-benchmarks directory.")
    parser.add_argument("--output-dir", default=Path(__file__).parent)
    args = parser.parse_args()

    from rdkit import Chem

    root = Path(args.cases).resolve().parent.parent
    data = Path(args.benchmarks) / "openfe_benchmarks/data/benchmark_systems"

    entries = []
    for case in json.loads(Path(args.cases).read_text()):
        (name_a, name_b), elements, mapping = openfe_mapping(root / case["case"])
        sdf = data / case["system"] / "ligands.sdf"
        ligands = {m.GetProp("_Name"): m for m in Chem.SDMolSupplier(str(sdf), removeHs=False)}
        for name, z in zip((name_a, name_b), elements):
            if [a.GetAtomicNum() for a in ligands[name].GetAtoms()] != z:
                raise ValueError(f"{case['case']}: atom order of {name} differs from {sdf}")
        entries.append(
            {
                "system": case["system"],
                "id": f"{name_a}~{name_b}",
                "ligand_a": name_a,
                "ligand_b": name_b,
                "mapping": {str(k): v for k, v in mapping.items()},
                "metadata": {"mode": case["mode"]},
                "source": case["case"],
            }
        )
        print(f"{case['case']}: {name_a}~{name_b}, {len(mapping)} mapped atoms")

    output_dir = Path(args.output_dir)
    (output_dir / "mappings_openfe.json").write_text(json.dumps(entries, indent=1))
    bss = [{**entry, "mapper": "biosimspace"} for entry in entries]
    (output_dir / "mappings_biosimspace.json").write_text(json.dumps(bss, indent=1))


if __name__ == "__main__":
    main()
