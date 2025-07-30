import os

import requests


def download_file(url, filepath, verbose=True):
    if verbose:
        print(f"Downloading {url} to {filepath}")
    response = requests.get(url)

    target_dir = os.path.dirname(filepath)
    if target_dir and not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Check if the request was successful
    if response.status_code == 200:
        with open(filepath, "wb") as file:
            file.write(response.content)
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

    return filepath
