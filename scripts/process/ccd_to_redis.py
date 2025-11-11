import argparse
import pickle

from rdkit import Chem
from redis import Redis


def load_pickled_dict_into_redis(
    path_in: str,
    host: str = "localhost",
    port: int = 7777,
    db: int = 0,
    batch: int = 10_000,
):
    """
    Load a pickled dict of {str -> any_serializable} into Redis.

    Parameters
    ----------
    path_in : str
        Path to the pickled dict file.
    host : str, optional
        Redis host, by default "localhost".
    port : int, optional
        Redis port, by default 7777.
    db : int, optional
        Redis database number, by default 0.
    batch : int, optional
        Number of keys to set per pipeline execution, by default 10_000.
    """
    # 1) read your dict
    with open(path_in, "rb") as f:
        data = pickle.load(f)  # {str -> any_serializable}

    # 2) connect; keep decode_responses=False so .get() returns bytes
    r = Redis(host=host, port=port, db=db, decode_responses=False)

    # 3) write as pickled bytes
    # we want all RDKit properties pickled
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    pipe = r.pipeline(transaction=False)
    i = 0
    for k, v in data.items():
        if not isinstance(k, str):
            raise TypeError(f"Key is not str: {k!r}")
        b = pickle.dumps(v, protocol=pickle.HIGHEST_PROTOCOL)
        pipe.set(k, b)
        i += 1
        if i % batch == 0:
            pipe.execute()
    pipe.execute()

    # 4) Snapshot to path_out
    r.bgsave()


if __name__ == "__main__":
    """
    Assumes Redis is running locally.

    ./redis-server \
      --save "" --appendonly no --port 7777 \
      --dir /abs/path/to/outdir \
      --dbfilename mydump.rdb
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--ccd_in", type=str)
    args = parser.parse_args()
    load_pickled_dict_into_redis(
        path_in=args.ccd_in,
    )
