import os

import httpx
from tqdm import tqdm


def download_file(url, filepath, verbose=True):
    target_dir = os.path.dirname(filepath)
    if target_dir and not os.path.exists(target_dir):
        os.makedirs(target_dir)

    with httpx.stream("GET", url, follow_redirects=True) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length", 0))
        if verbose:
            print(f"Downloading {url} to {filepath}")
            with open(filepath, "wb") as file, tqdm(
                total=total or None,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=os.path.basename(filepath),
            ) as pbar:
                for chunk in response.iter_bytes(chunk_size=8192):
                    file.write(chunk)
                    pbar.update(len(chunk))
        else:
            with open(filepath, "wb") as file:
                for chunk in response.iter_bytes(chunk_size=8192):
                    file.write(chunk)

    return filepath
