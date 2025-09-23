import argparse
import concurrent.futures
import subprocess
from pathlib import Path

from tqdm import tqdm

OST_COMPARE_STRUCTURE = r"""
#!/bin/bash
# https://openstructure.org/docs/2.7/actions/#ost-compare-structures

IMAGE_NAME=openstructure-0.2.8

command="compare-structures \
-m {model_file} \
-r {reference_file} \
--fault-tolerant \
--min-pep-length 4 \
--min-nuc-length 4 \
-o {output_path} \
--lddt --bb-lddt --qs-score --dockq \
--ics --ips --rigid-scores --patch-scores --tm-score"

sudo docker run -u $(id -u):$(id -g) --rm --volume {mount}:{mount} $IMAGE_NAME $command
"""


OST_COMPARE_LIGAND = r"""
#!/bin/bash
# https://openstructure.org/docs/2.7/actions/#ost-compare-structures

IMAGE_NAME=openstructure-0.2.8

command="compare-ligand-structures \
-m {model_file} \
-r {reference_file} \
--fault-tolerant \
--lddt-pli --rmsd \
--substructure-match \
-o {output_path}"

sudo docker run -u $(id -u):$(id -g) --rm --volume {mount}:{mount} $IMAGE_NAME $command
"""


def evaluate_structure(
    name: str,
    pred: Path,
    reference: Path,
    outdir: str,
    mount: str,
    executable: str = "/bin/bash",
) -> None:
    """Evaluate the structure."""
    # Evaluate polymer metrics
    out_path = Path(outdir) / f"{name}.json"

    if out_path.exists():
        print(  # noqa: T201
            f"Skipping recomputation of {name} as protein json file already exists"
        )
    else:
        subprocess.run(
            OST_COMPARE_STRUCTURE.format(
                model_file=str(pred),
                reference_file=str(reference),
                output_path=str(out_path),
                mount=mount,
            ),
            shell=True,  # noqa: S602
            check=False,
            executable=executable,
            capture_output=True,
        )

    # Evaluate ligand metrics
    out_path = Path(outdir) / f"{name}_ligand.json"
    if out_path.exists():
        print(f"Skipping recomputation of {name} as ligand json file already exists")  # noqa: T201
    else:
        subprocess.run(
            OST_COMPARE_LIGAND.format(
                model_file=str(pred),
                reference_file=str(reference),
                output_path=str(out_path),
                mount=mount,
            ),
            shell=True,  # noqa: S602
            check=False,
            executable=executable,
            capture_output=True,
        )


def main(args):
    # Aggregate the predictions and references
    files = list(args.data.iterdir())
    names = {f.stem.lower(): f for f in files}

    # Create the output directory
    args.outdir.mkdir(parents=True, exist_ok=True)

    first_item = True
    with concurrent.futures.ThreadPoolExecutor(args.max_workers) as executor:
        futures = []
        for name, folder in names.items():
            for model_id in range(5):
                # Split the input data
                if args.format == "af3":
                    pred_path = folder / f"seed-1_sample-{model_id}" / "model.cif"
                elif args.format == "chai":
                    pred_path = folder / f"pred.model_idx_{model_id}.cif"
                elif args.format == "boltz":
                    name_file = (
                        f"{name[0].upper()}{name[1:]}"
                        if args.testset == "casp"
                        else name.lower()
                    )
                    pred_path = folder / f"{name_file}_model_{model_id}.cif"

                if args.testset == "casp":
                    ref_path = args.pdb / f"{name[0].upper()}{name[1:]}.cif"
                elif args.testset == "test":
                    ref_path = args.pdb / f"{name.lower()}.cif.gz"

                if first_item:
                    # Evaluate the first item in the first prediction
                    # Ensures that the docker image is downloaded
                    evaluate_structure(
                        name=f"{name}_model_{model_id}",
                        pred=str(pred_path),
                        reference=str(ref_path),
                        outdir=str(args.outdir),
                        mount=args.mount,
                        executable=args.executable,
                    )
                    first_item = False
                else:
                    future = executor.submit(
                        evaluate_structure,
                        name=f"{name}_model_{model_id}",
                        pred=str(pred_path),
                        reference=str(ref_path),
                        outdir=str(args.outdir),
                        mount=args.mount,
                        executable=args.executable,
                    )
                    futures.append(future)

        # Wait for all tasks to complete
        with tqdm(total=len(futures)) as pbar:
            for _ in concurrent.futures.as_completed(futures):
                pbar.update(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("pdb", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--format", type=str, default="af3")
    parser.add_argument("--testset", type=str, default="casp")
    parser.add_argument("--mount", type=str)
    parser.add_argument("--executable", type=str, default="/bin/bash")
    parser.add_argument("--max-workers", type=int, default=32)
    args = parser.parse_args()
    main(args)
