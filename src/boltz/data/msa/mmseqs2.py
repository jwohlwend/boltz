# From https://github.com/sokrypton/ColabFold/blob/main/colabfold/colabfold.py

import logging
import os
import random
import tarfile
import time
from typing import Optional, Union, Dict

import requests
from requests.auth import HTTPBasicAuth
from tqdm import tqdm

logger = logging.getLogger(__name__)

TQDM_BAR_FORMAT = (
    "{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]"
)


def run_mmseqs2(  # noqa: PLR0912, D103, C901, PLR0915
    x: Union[str, list[str]],
    prefix: str = "tmp",
    use_env: bool = True,
    use_filter: bool = True,
    use_pairing: bool = False,
    pairing_strategy: str = "greedy",
    host_url: str = "https://api.colabfold.com",
    msa_server_username: Optional[str] = None,
    msa_server_password: Optional[str] = None,
    auth_headers: Optional[Dict[str, str]] = None,
) -> tuple[list[str], list[str]]:
    submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

    # Validate mutually exclusive authentication methods
    has_basic_auth = msa_server_username and msa_server_password
    has_header_auth = auth_headers is not None
    if has_basic_auth and (has_header_auth or auth_headers):
        raise ValueError(
            "Cannot use both basic authentication (username/password) and header/API key authentication. "
            "Please use only one authentication method."
        )

    # Set header agent as boltz
    headers = {}
    headers["User-Agent"] = "boltz"

    # Set up authentication
    auth = None
    if has_basic_auth:
        auth = HTTPBasicAuth(msa_server_username, msa_server_password)
        logger.debug(f"MMSeqs2 server authentication: using basic auth for user '{msa_server_username}'")
    elif has_header_auth:
        headers.update(auth_headers)
        logger.debug("MMSeqs2 server authentication: using header-based authentication")
    else:
        logger.debug("MMSeqs2 server authentication: no credentials provided")
    
    logger.debug(f"Connecting to MMSeqs2 server at: {host_url}")
    logger.debug(f"Using endpoint: {submission_endpoint}")
    logger.debug(f"Pairing strategy: {pairing_strategy}")
    logger.debug(f"Use environment databases: {use_env}")
    logger.debug(f"Use filtering: {use_filter}")

    def submit(seqs, mode, N=101):
        n, query = N, ""
        for seq in seqs:
            query += f">{n}\n{seq}\n"
            n += 1

        error_count = 0
        while True:
            try:
                # https://requests.readthedocs.io/en/latest/user/advanced/#advanced
                # "good practice to set connect timeouts to slightly larger than a multiple of 3"
                logger.debug(f"Submitting MSA request to {host_url}/{submission_endpoint}")
                res = requests.post(
                    f"{host_url}/{submission_endpoint}",
                    data={"q": query, "mode": mode},
                    timeout=6.02,
                    headers=headers,
                    auth=auth,
                )
                logger.debug(f"MSA submission response status: {res.status_code}")
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                if error_count > 5:
                    raise Exception(
                        "Too many failed attempts for the MSA generation request."
                    )
                time.sleep(5)
            else:
                break

        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def status(ID):
        error_count = 0
        while True:
            try:
                logger.debug(f"Checking MSA job status for ID: {ID}")
                res = requests.get(
                    f"{host_url}/ticket/{ID}", timeout=6.02, headers=headers, auth=auth
                )
                logger.debug(f"MSA status check response status: {res.status_code}")
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                if error_count > 5:
                    raise Exception(
                        "Too many failed attempts for the MSA generation request."
                    )
                time.sleep(5)
            else:
                break
        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def download(ID, path):
        error_count = 0
        while True:
            try:
                logger.debug(f"Downloading MSA results for ID: {ID}")
                res = requests.get(
                    f"{host_url}/result/download/{ID}", timeout=6.02, headers=headers, auth=auth
                )
                logger.debug(f"MSA download response status: {res.status_code}")
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                if error_count > 5:
                    raise Exception(
                        "Too many failed attempts for the MSA generation request."
                    )
                time.sleep(5)
            else:
                break
        with open(path, "wb") as out:
            out.write(res.content)

    # process input x
    seqs = [x] if isinstance(x, str) else x

    # setup mode
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    if use_pairing:
        mode = ""
        # greedy is default, complete was the previous behavior
        if pairing_strategy == "greedy":
            mode = "pairgreedy"
        elif pairing_strategy == "complete":
            mode = "paircomplete"
        if use_env:
            mode = mode + "-env"

    # define path
    path = f"{prefix}_{mode}"
    if not os.path.isdir(path):
        os.mkdir(path)

    # call mmseqs2 api
    tar_gz_file = f"{path}/out.tar.gz"
    N, REDO = 101, True

    # deduplicate and keep track of order
    seqs_unique = []
    # TODO this might be slow for large sets
    [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    # lets do it!
    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 150 * len(seqs_unique)
        with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = submit(seqs_unique, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    sleep_time = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
                    # resubmit
                    time.sleep(sleep_time)
                    out = submit(seqs_unique, mode, N)

                if out["status"] == "ERROR":
                    msg = (
                        "MMseqs2 API is giving errors. Please confirm your "
                        " input is a valid protein sequence. If error persists, "
                        "please try again an hour later."
                    )
                    raise Exception(msg)

                if out["status"] == "MAINTENANCE":
                    msg = (
                        "MMseqs2 API is undergoing maintenance. "
                        "Please try again in a few minutes."
                    )
                    raise Exception(msg)

                # wait for job to finish
                ID, TIME = out["id"], 0
                logger.debug(f"MSA job submitted successfully with ID: {ID}")
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                        pbar.update(n=t)

                if out["status"] == "COMPLETE":
                    logger.debug(f"MSA job completed successfully for ID: {ID}")
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE - TIME))
                    REDO = False

                if out["status"] == "ERROR":
                    REDO = False
                    msg = (
                        "MMseqs2 API is giving errors. Please confirm your "
                        " input is a valid protein sequence. If error persists, "
                        "please try again an hour later."
                    )
                    raise Exception(msg)

            # Download results
            download(ID, tar_gz_file)

    # prep list of a3m files
    if use_pairing:
        a3m_files = [f"{path}/pair.a3m"]
    else:
        a3m_files = [f"{path}/uniref.a3m"]
        if use_env:
            a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

    # extract a3m files
    if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(path)

    # gather a3m lines
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        for line in open(a3m_file, "r"):
            if len(line) > 0:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    M = int(line[1:].rstrip())
                    update_M = False
                    if M not in a3m_lines:
                        a3m_lines[M] = []
                a3m_lines[M].append(line)

    a3m_lines = ["".join(a3m_lines[n]) for n in Ms]
    return a3m_lines
