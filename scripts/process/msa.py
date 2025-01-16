import argparse
import multiprocessing
from dataclasses import asdict
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np
from p_tqdm import p_umap
from redis import Redis
from tqdm import tqdm

from boltz.data.parse.a3m import parse_a3m


class Resource:
    """A shared resource for processing."""

    def __init__(self, host: str, port: int) -> None:
        """Initialize the redis database."""
        self._redis = Redis(host=host, port=port)

    def get(self, key: str) -> Any:  # noqa: ANN401
        """Get an item from the Redis database."""
        return self._redis.get(key)

    def __getitem__(self, key: str) -> Any:  # noqa: ANN401
        """Get an item from the resource."""
        out = self.get(key)
        if out is None:
            raise KeyError(key)
        return out


def process_msa(
    path: Path,
    outdir: str,
    max_seqs: int,
    resource: Resource,
) -> None:
    """Run processing in a worker thread."""
    outdir = Path(outdir)
    out_path = outdir / f"{path.stem}.npz"
    if not out_path.exists():
        msa = parse_a3m(path, resource, max_seqs)
        np.savez_compressed(out_path, **asdict(msa))


def process(args) -> None:
    """Run the data processing task."""
    # Create output directory
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load the resource
    resource = Resource(host=args.redis_host, port=args.redis_port)

    # Get data points
    print("Fetching data...")
    data = list(args.msadir.rglob("*.a3m*"))
    print(f"Found {len(data)} MSA's.")

    # Check if we can run in parallel
    max_processes = multiprocessing.cpu_count()
    num_processes = max(1, min(args.num_processes, max_processes, len(data)))
    parallel = num_processes > 1

    # Run processing
    if parallel:
        # Create processing function
        fn = partial(
            process_msa,
            outdir=args.outdir,
            max_seqs=args.max_seqs,
            resource=resource,
        )

        # Run in parallel
        p_umap(fn, data, num_cpus=num_processes)

    else:
        # Run in serial
        for path in tqdm(data):
            process_msa(
                path,
                outdir=args.outdir,
                max_seqs=args.max_seqs,
                resource=resource,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MSA data.")
    parser.add_argument(
        "--msadir",
        type=Path,
        required=True,
        help="The MSA data directory.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default="data",
        help="The output directory.",
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        default=multiprocessing.cpu_count(),
        help="The number of processes.",
    )
    parser.add_argument(
        "--redis-host",
        type=str,
        default="localhost",
        help="The Redis host.",
    )
    parser.add_argument(
        "--redis-port",
        type=int,
        default=7777,
        help="The Redis port.",
    )
    parser.add_argument(
        "--max-seqs",
        type=int,
        default=16384,
        help="The maximum number of sequences.",
    )
    args = parser.parse_args()
    process(args)
