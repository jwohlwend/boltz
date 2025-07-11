from pyrosetta import *
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task import TaskFactory

from pathlib import Path
import json
import pandas as pd
from tqdm.auto import tqdm

import sys
from contextlib import redirect_stdout, redirect_stderr


def ensure_pyrosetta_initialized():
    """
    Ensure that PyRosetta is initialized.
    """
    if not hasattr(ensure_pyrosetta_initialized, "_initialized"):
        init(silent=False, set_logging_handler="logging")
        ensure_pyrosetta_initialized._initialized = True


def get_rosetta_score(input_path):
    """
    Computes the Rosetta score.
    Args:
        input_path (str): Path to the CIF/PDB file.
    Returns:
        float: Rosetta energy score.
    """
    ensure_pyrosetta_initialized()

    try:
        pose = pose_from_file(str(input_path))
        # Get the ref2015 score function
        scorefxn = get_score_function()
        return scorefxn(pose)
    except:
        raise RuntimeError(f"ERROR in {input_path}!")


def repack_sidechains(input_path, out_path, **kwargs):
    """
        Performs Rosetta side-chain repacking and returns the final score
    Args:
        input_path (str): Path to the CIF/PDB file.
        out_path (str): Path to the output PDB file.
    Returns:
        float: repacked Rosetta score.
    """
    ensure_pyrosetta_initialized()

    pose = pose_from_file(str(input_path))
    scorefxn = get_score_function()
    tf = TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
    packer = PackRotamersMover()
    packer.task_factory(tf)
    packer.score_function(scorefxn)
    packer.apply(pose)
    pose.dump_pdb(str(out_path))
    return scorefxn(pose)


def fastrelax(input_path, out_path, constrain_relax_to_start_coords=True, **kwargs):
    """
    Performs Rosetta FastRelax and returns the final score
    Args:
        input_path (str): Path to the CIF/PDB file.
        out_path (str): Path to the output PDB file.
    Returns:
        float: fastrelaxed Rosetta score.
    """
    ensure_pyrosetta_initialized()

    pose = pose_from_file(str(input_path))
    scorefxn = get_score_function()
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(scorefxn)
    fast_relax.constrain_relax_to_start_coords(constrain_relax_to_start_coords)
    fast_relax.apply(pose)
    pose.dump_pdb(str(out_path))
    return scorefxn(pose)


def relax(
    input_path,
    output_dir=None,
    override=False,
    save_logs=True,
    save_energies=True,
    **kwargs,
):
    input_path = Path(input_path)
    if output_dir is None:
        output_dir = input_path.parent
    stdout = (
        str(output_dir / f"rosetta-relax_{input_path.stem}.stdout")
        if save_logs
        else "/dev/null"
    )
    stderr = (
        str(output_dir / f"rosetta-relax_{input_path.stem}.stderr")
        if save_logs
        else "/dev/null"
    )

    with redirect_stdout(open(stdout, "w")):
        with redirect_stderr(open(stderr, "w")):
            init_energy = get_rosetta_score(input_path)
            ret = dict(init_energy=init_energy, input_path=str(input_path))

            repacked_path = output_dir / f"repacked_{input_path.stem}.pdb"
            repacked_energy = (
                repack_sidechains(input_path, repacked_path, **kwargs)
                if not repacked_path.exists() or override
                else get_rosetta_score(repacked_path)
            )
            ret.update(
                dict(repacked_energy=repacked_energy, repacked_path=str(repacked_path))
            )

            relax_path = output_dir / f"fastrelaxed_{input_path.stem}.pdb"
            fastrelaxed_energy = (
                fastrelax(repacked_path, relax_path, **kwargs)
                if not relax_path.exists() or override
                else get_rosetta_score(relax_path)
            )
            ret.update(
                dict(
                    fastrelaxed_energy=fastrelaxed_energy,
                    fastrelaxed_path=str(relax_path),
                )
            )

            if save_energies:
                json_path = output_dir / f"energies_{input_path.stem}.json"
                with open(json_path, "w") as f:
                    json.dump(
                        ret, f, indent=4
                    )  # `indent=4` makes the file human-readable

            sys.stdout.flush()
            sys.stderr.flush()

    return ret


def parallel_relax(
    input_paths,
    output_dir=None,
    override=False,
    save_logs=False,
    save_energies=True,
    cores=8,
    **kwargs,
):
    import multiprocessing
    from functools import partial

    if not all([Path(p).suffix.lower() in [".cif", ".mmcif"] for p in input_paths]):
        print(
            "WARNING: If your structure contains a ligand, use CCD ligands in Boltz input and choose CIF as the output format. "
            "This will ensure that the atom/residue naming is compatible with Rosetta."
        )
    with multiprocessing.get_context("spawn").Pool(cores) as ex:
        ret = pd.DataFrame(
            tqdm(
                ex.imap_unordered(
                    partial(
                        relax,
                        output_dir=output_dir,
                        override=override,
                        save_logs=save_logs,
                        save_energies=save_energies,
                        **kwargs,
                    ),
                    input_paths,
                ),
                total=len(input_paths),
                desc="Rosetta Relaxation",
            )
        )
    ret["name"] = [Path(p).parent.name for p in ret.input_path]
    return ret
