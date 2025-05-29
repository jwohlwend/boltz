"""Compute conformers and symmetries for all the CCD molecules."""

import argparse
import multiprocessing
import pickle
import sys
from functools import partial
from pathlib import Path

import pandas as pd
import rdkit
from p_tqdm import p_uimap
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core.component import ConformerType
from rdkit import rdBase
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Conformer, Mol
from tqdm import tqdm


def load_molecules(components: str) -> list[Mol]:
    """Load the CCD components file.

    Parameters
    ----------
    components : str
        Path to the CCD components file.

    Returns
    -------
    list[Mol]

    """
    components: dict[str, ccd_reader.CCDReaderResult]
    components = ccd_reader.read_pdb_components_file(components)

    mols = []
    for name, component in components.items():
        mol = component.component.mol
        mol.SetProp("PDB_NAME", name)
        mols.append(mol)

    return mols


def compute_3d(mol: Mol, version: str = "v3") -> bool:
    """Generate 3D coordinates using EKTDG method.

    Taken from `pdbeccdutils.core.component.Component`.

    Parameters
    ----------
    mol: Mol
        The RDKit molecule to process
    version: str, optional
        The ETKDG version, defaults ot v3

    Returns
    -------
    bool
        Whether computation was successful.

    """
    if version == "v3":
        options = rdkit.Chem.AllChem.ETKDGv3()
    elif version == "v2":
        options = rdkit.Chem.AllChem.ETKDGv2()
    else:
        options = rdkit.Chem.AllChem.ETKDGv2()

    options.clearConfs = False
    conf_id = -1

    try:
        conf_id = rdkit.Chem.AllChem.EmbedMolecule(mol, options)
        rdkit.Chem.AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=1000)

    except RuntimeError:
        pass  # Force field issue here
    except ValueError:
        pass  # sanitization issue here

    if conf_id != -1:
        conformer = mol.GetConformer(conf_id)
        conformer.SetProp("name", ConformerType.Computed.name)
        conformer.SetProp("coord_generation", f"ETKDG{version}")

        return True

    return False


def get_conformer(mol: Mol, c_type: ConformerType) -> Conformer:
    """Retrieve an rdkit object for a deemed conformer.

    Taken from `pdbeccdutils.core.component.Component`.

    Parameters
    ----------
    mol: Mol
        The molecule to process.
    c_type: ConformerType
        The conformer type to extract.

    Returns
    -------
    Conformer
        The desired conformer, if any.

    Raises
    ------
    ValueError
        If there are no conformers of the given tyoe.

    """
    for c in mol.GetConformers():
        try:
            if c.GetProp("name") == c_type.name:
                return c
        except KeyError:  # noqa: PERF203
            pass

    msg = f"Conformer {c_type.name} does not exist."
    raise ValueError(msg)


def compute_symmetries(mol: Mol) -> list[list[int]]:
    """Compute the symmetries of a molecule.

    Parameters
    ----------
    mol : Mol
        The molecule to process

    Returns
    -------
    list[list[int]]
        The symmetries as a list of index permutations

    """
    mol = AllChem.RemoveHs(mol)
    idx_map = {}
    atom_idx = 0
    for i, atom in enumerate(mol.GetAtoms()):
        # Skip if leaving atoms
        if int(atom.GetProp("leaving_atom")):
            continue
        idx_map[i] = atom_idx
        atom_idx += 1

    # Calculate self permutations
    permutations = []
    raw_permutations = mol.GetSubstructMatches(mol, uniquify=False)
    for raw_permutation in raw_permutations:
        # Filter out permutations with leaving atoms
        try:
            if {raw_permutation[idx] for idx in idx_map} == set(idx_map.keys()):
                permutation = [
                    idx_map[idx] for idx in raw_permutation if idx in idx_map
                ]
                permutations.append(permutation)
        except Exception:  # noqa: S110, PERF203, BLE001
            pass
    serialized_permutations = pickle.dumps(permutations)
    mol.SetProp("symmetries", serialized_permutations.hex())
    return permutations


def process(mol: Mol, output: str) -> tuple[str, str]:
    """Process a CCD component.

    Parameters
    ----------
    mol : Mol
        The molecule to process
    output : str
        The directory to save the molecules

    Returns
    -------
    str
        The name of the component
    str
        The result of the conformer generation

    """
    # Get name
    name = mol.GetProp("PDB_NAME")

    # Check if single atom
    if mol.GetNumAtoms() == 1:
        result = "single"
    else:
        # Get the 3D conformer
        try:
            # Try to generate a 3D conformer with RDKit
            success = compute_3d(mol, version="v3")
            if success:
                _ = get_conformer(mol, ConformerType.Computed)
                result = "computed"

            # Otherwise, default to the ideal coordinates
            else:
                _ = get_conformer(mol, ConformerType.Ideal)
                result = "ideal"
        except ValueError:
            result = "failed"

    # Dump the molecule
    path = Path(output) / f"{name}.pkl"
    with path.open("wb") as f:
        pickle.dump(mol, f)

    # Output the results
    return name, result


def main(args: argparse.Namespace) -> None:
    """Process conformers."""
    # Set property saving
    rdkit.Chem.SetDefaultPickleProperties(rdkit.Chem.PropertyPickleOptions.AllProps)

    # Load components
    print("Loading components")  # noqa: T201
    molecules = load_molecules(args.components)

    # Reset stdout and stderr, as pdbccdutils messes with them
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    # Disable rdkit warnings
    blocker = rdBase.BlockLogs()  # noqa: F841

    # Setup processing function
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mol_output = outdir / "mols"
    mol_output.mkdir(parents=True, exist_ok=True)
    process_fn = partial(process, output=str(mol_output))

    # Process the files in parallel
    print("Processing components")  # noqa: T201
    metadata = []

    # Check if we can run in parallel
    max_processes = multiprocessing.cpu_count()
    num_processes = max(1, min(args.num_processes, max_processes, len(molecules)))
    parallel = num_processes > 1

    if parallel:
        for name, result in p_uimap(
            process_fn,
            molecules,
            num_cpus=num_processes,
        ):
            metadata.append({"name": name, "result": result})
    else:
        for mol in tqdm(molecules):
            name, result = process_fn(mol)
            metadata.append({"name": name, "result": result})

    # Load and group outputs
    molecules = {}
    for item in metadata:
        if item["result"] == "failed":
            continue

        # Load the mol file
        path = mol_output / f"{item['name']}.pkl"
        with path.open("rb") as f:
            mol = pickle.load(f)  # noqa: S301
            molecules[item["name"]] = mol

    # Dump metadata
    path = outdir / "results.csv"
    metadata = pd.DataFrame(metadata)
    metadata.to_csv(path)

    # Dump the components
    path = outdir / "ccd.pkl"
    with path.open("wb") as f:
        pickle.dump(molecules, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--components", type=str)
    parser.add_argument("--outdir", type=str)
    parser.add_argument(
        "--num_processes",
        type=int,
        default=multiprocessing.cpu_count(),
    )
    args = parser.parse_args()
    main(args)
