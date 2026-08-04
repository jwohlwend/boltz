from __future__ import annotations

import os
import pickle
import platform
import re
import tarfile
import warnings
from dataclasses import asdict, dataclass, replace
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import TYPE_CHECKING, List, Literal, Optional, Union

import click

from importlib.metadata import version as _pkg_version

try:
    from pytorch_lightning.utilities import rank_zero_only
except ImportError:

    def rank_zero_only(fn):  # type: ignore[misc]
        return fn


if TYPE_CHECKING:
    from boltz.data.types import Manifest

CCD_URL = "https://huggingface.co/boltz-community/boltz-1/resolve/main/ccd.pkl"
MOL_URL = "https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar"

BOLTZ1_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz1_conf.ckpt",
    "https://huggingface.co/boltz-community/boltz-1/resolve/main/boltz1_conf.ckpt",
]

BOLTZ2_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz2_conf.ckpt",
    "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt",
]

BOLTZ2_AFFINITY_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz2_aff.ckpt",
    "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_aff.ckpt",
]


def _resolved_path_for_message(path: Path) -> Path:
    """Return a best-effort resolved path for user-facing error messages."""
    try:
        return path.expanduser().resolve()
    except OSError:
        return path.expanduser()


def _checkpoint_needs_download(path: Path, description: str) -> bool:
    """Return True when a cached checkpoint should be downloaded again."""
    resolved = _resolved_path_for_message(path)

    if not path.exists():
        return True
    if not path.is_file():
        msg = f"{description} checkpoint path is not a file: {resolved}"
        raise RuntimeError(msg)

    try:
        size = path.stat().st_size
    except OSError as e:
        msg = f"Could not read {description} checkpoint at {resolved}: {e}"
        raise RuntimeError(msg) from e

    if size == 0:
        click.echo(f"Found empty {description} checkpoint at {resolved}; re-downloading.")
        return True

    return False


def _validate_checkpoint_for_load(path: Path, description: str) -> None:
    """Raise a clear CLI error for missing or empty checkpoint files."""
    resolved = _resolved_path_for_message(path)

    if not path.exists():
        msg = (
            f"{description} checkpoint not found at {resolved}. "
            "Delete any partial file and retry so boltz can re-download it, "
            "or pass a valid checkpoint path."
        )
        raise click.ClickException(msg)
    if not path.is_file():
        msg = f"{description} checkpoint path is not a file: {resolved}"
        raise click.ClickException(msg)

    try:
        size = path.stat().st_size
    except OSError as e:
        msg = f"Could not read {description} checkpoint at {resolved}: {e}"
        raise click.ClickException(msg) from e

    if size == 0:
        msg = (
            f"{description} checkpoint is empty: {resolved}. "
            "Delete the file and retry so boltz can re-download it, "
            "or pass a valid checkpoint path."
        )
        raise click.ClickException(msg)


def _available_cpu_count() -> int:
    """Return the number of CPUs available to this process.

    Respects cgroup/taskset limits on Linux via os.sched_getaffinity.
    Falls back to os.cpu_count on macOS or on error.
    """
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def _resolve_affinity_max_parallel_samples(
    max_parallel_samples: Optional[int],
    diffusion_samples_affinity: int,
) -> int:
    """Resolve max_parallel_samples for the affinity diffusion path.

    Upstream hardcoded this to 1, which was silently a no-op under upstream's
    buggy chunking math (``chunk(multiplicity % max + 1)`` collapsed to
    ``chunk(1)`` and batched all samples into one pass). Our diffusion fix
    made the hardcode actually take effect, which forced N sequential
    single-sample forward passes per record — an avoidable ~Nx slowdown.

    We now honor the user's ``--max_parallel_samples``, capped at the number
    of affinity diffusion samples so the default cleanly batches all samples
    in one pass without inflating the work above what diffusion will run.

    ``max_parallel_samples=None`` is treated the same way the diffusion code
    treats it ("use multiplicity") — i.e. batch all affinity samples in one
    pass. The Click CLI default is 5, so this branch is only reached by
    programmatic callers of ``predict.callback(...)`` that pass ``None``.
    """
    if max_parallel_samples is None:
        return diffusion_samples_affinity
    return min(max_parallel_samples, diffusion_samples_affinity)


def _resolve_checkpoint_load_map_location(accelerator: str) -> str:
    """Return the checkpoint restore device for the requested accelerator.

    Some newer PyTorch MPS builds advertise Metal support but do not expose
    ``torch.mps.current_device()``, which breaks checkpoint deserialization when
    ``map_location='mps'`` is passed into ``torch.load``/Lightning. In that case,
    we fall back to CPU for loading and let the Trainer move the model to MPS
    after it has been created.
    """
    if accelerator == "cpu":
        return "cpu"
    if accelerator == "mps":
        try:
            import torch
        except ImportError:
            return "cpu"

        mps_available = getattr(torch.backends, "mps", None)
        if mps_available is None or not mps_available.is_available():
            return "cpu"

        if not hasattr(torch.mps, "current_device"):
            return "cpu"
        return "mps"
    return "cuda"


def _ensure_torch_mps_current_device_compat() -> None:
    """Provide a minimal ``torch.mps.current_device`` shim when missing.

    Some nightly/dev builds report MPS availability but omit the
    ``current_device`` attribute entirely. That breaks Torch's deserializer
    when it probes the active Metal device while restoring checkpoint tensors.
    """
    try:
        import torch
    except ImportError:
        return

    if not hasattr(torch.mps, "current_device"):
        torch.mps.current_device = lambda: 0  # type: ignore[attr-defined]


def _normalize_accelerator(accelerator: str) -> str:
    """Use MPS for Apple Silicon when the user requests the generic GPU path."""
    try:
        import torch
    except ImportError:
        return accelerator

    if accelerator == "gpu" and getattr(torch.backends, "mps", None) is not None:
        if torch.backends.mps.is_available():
            return "mps"
    return accelerator


@dataclass
class BoltzProcessedInput:
    """Processed input data."""

    manifest: Manifest
    targets_dir: Path
    msa_dir: Path
    constraints_dir: Optional[Path] = None
    template_dir: Optional[Path] = None
    extra_mols_dir: Optional[Path] = None


@dataclass
class PairformerArgs:
    """Pairformer arguments."""

    num_blocks: int = 48
    num_heads: int = 16
    dropout: float = 0.0
    activation_checkpointing: bool = False
    v2: bool = False


@dataclass
class PairformerArgsV2:
    """Pairformer arguments."""

    num_blocks: int = 64
    num_heads: int = 16
    dropout: float = 0.0
    activation_checkpointing: bool = False


@dataclass
class MSAModuleArgs:
    """MSA module arguments."""

    msa_s: int = 64
    msa_blocks: int = 4
    msa_dropout: float = 0.0
    z_dropout: float = 0.0
    use_paired_feature: bool = True
    pairwise_head_width: int = 32
    pairwise_num_heads: int = 4
    activation_checkpointing: bool = False
    subsample_msa: bool = False
    num_subsampled_msa: int = 1024


@dataclass
class BoltzDiffusionParams:
    """Diffusion process parameters."""

    gamma_0: float = 0.605
    gamma_min: float = 1.107
    noise_scale: float = 0.901
    rho: float = 8
    step_scale: float = 1.638
    sigma_min: float = 0.0004
    sigma_max: float = 160.0
    sigma_data: float = 16.0
    P_mean: float = -1.2
    P_std: float = 1.5
    coordinate_augmentation: bool = True
    alignment_reverse_diff: bool = True
    synchronize_sigmas: bool = True
    use_inference_model_cache: bool = True


@dataclass
class Boltz2DiffusionParams:
    """Diffusion process parameters."""

    gamma_0: float = 0.8
    gamma_min: float = 1.0
    noise_scale: float = 1.003
    rho: float = 7
    step_scale: float = 1.5
    sigma_min: float = 0.0001
    sigma_max: float = 160.0
    sigma_data: float = 16.0
    P_mean: float = -1.2
    P_std: float = 1.5
    coordinate_augmentation: bool = True
    alignment_reverse_diff: bool = True
    synchronize_sigmas: bool = True


@dataclass
class BoltzSteeringParams:
    """Steering parameters."""

    fk_steering: bool = False
    num_particles: int = 3
    fk_lambda: float = 4.0
    fk_resampling_interval: int = 3
    physical_guidance_update: bool = False
    contact_guidance_update: bool = True
    num_gd_steps: int = 20


def _download(url: str, dest: str) -> None:
    """Download a file with proper redirect handling.

    Uses requests instead of urllib to handle HuggingFace redirects
    and to use certifi SSL certificates (avoids macOS SSL errors).
    """
    import requests

    response = requests.get(url, stream=True, allow_redirects=True, timeout=300)
    response.raise_for_status()
    with open(dest, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)


@rank_zero_only
def download_boltz1(cache: Path) -> None:
    """Download all the required data.

    Parameters
    ----------
    cache : Path
        The cache directory.

    """
    # Download CCD
    ccd = cache / "ccd.pkl"
    if not ccd.exists():
        click.echo(
            f"Downloading the CCD dictionary to {ccd}. You may "
            "change the cache directory with the --cache flag."
        )
        _download(CCD_URL, str(ccd))
    # Download model
    model = cache / "boltz1_conf.ckpt"
    if _checkpoint_needs_download(model, "Boltz-1 weights"):
        click.echo(
            f"Downloading the model weights to {model}. You may "
            "change the cache directory with the --cache flag."
        )
        for i, url in enumerate(BOLTZ1_URL_WITH_FALLBACK):
            try:
                _download(url, str(model))
                _validate_checkpoint_for_load(model, "Boltz-1 weights")
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ1_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                continue


@rank_zero_only
def download_boltz2(cache: Path) -> None:
    """Download all the required data.

    Parameters
    ----------
    cache : Path
        The cache directory.

    """
    # Download CCD
    mols = cache / "mols"
    tar_mols = cache / "mols.tar"
    if not mols.exists():
        if not tar_mols.exists():
            click.echo(
                f"Downloading the CCD data to {tar_mols}. "
                "This may take a bit of time. You may change the cache directory "
                "with the --cache flag."
            )
            _download(MOL_URL, str(tar_mols))
        click.echo(
            f"Extracting the CCD data to {mols}. "
            "This may take a bit of time. You may change the cache directory "
            "with the --cache flag."
        )
        with tarfile.open(str(tar_mols), "r") as tar:
            tar.extractall(cache)  # noqa: S202

    # Download model
    model = cache / "boltz2_conf.ckpt"
    if _checkpoint_needs_download(model, "Boltz-2 weights"):
        click.echo(
            f"Downloading the Boltz-2 weights to {model}. You may "
            "change the cache directory with the --cache flag."
        )
        for i, url in enumerate(BOLTZ2_URL_WITH_FALLBACK):
            try:
                _download(url, str(model))
                _validate_checkpoint_for_load(model, "Boltz-2 weights")
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ2_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                continue

    # Download affinity model
    affinity_model = cache / "boltz2_aff.ckpt"
    if _checkpoint_needs_download(affinity_model, "Boltz-2 affinity weights"):
        click.echo(
            f"Downloading the Boltz-2 affinity weights to {affinity_model}. You may "
            "change the cache directory with the --cache flag."
        )
        for i, url in enumerate(BOLTZ2_AFFINITY_URL_WITH_FALLBACK):
            try:
                _download(url, str(affinity_model))
                _validate_checkpoint_for_load(
                    affinity_model, "Boltz-2 affinity weights"
                )
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ2_AFFINITY_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                continue


def get_cache_path() -> str:
    """Determine the cache path, prioritising the BOLTZ_CACHE environment variable.

    Returns
    -------
    str: Path
        Path to use for boltz cache location.

    """
    env_cache = os.environ.get("BOLTZ_CACHE")
    if env_cache:
        resolved_cache = Path(env_cache).expanduser().resolve()
        if not resolved_cache.is_absolute():
            msg = f"BOLTZ_CACHE must be an absolute path, got: {env_cache}"
            raise ValueError(msg)
        return str(resolved_cache)

    return str(Path("~/.boltz").expanduser())


def check_inputs(data: Path, skip_bad_inputs: bool = False) -> list[Path]:
    """Check the input data and output directory.

    Parameters
    ----------
    data : Path
        The input data.
    skip_bad_inputs : bool
        If True, filter out invalid entries and warn instead of aborting.

    Returns
    -------
    list[Path]
        The list of input data.

    """
    click.echo("Checking input data.")

    # Check if data is a directory
    if data.is_dir():
        data: list[Path] = list(data.glob("*"))

        # Filter out non .fasta or .yaml files
        valid = []
        for d in data:
            if d.is_dir():
                msg = f"Found directory {d} instead of .fasta or .yaml."
            elif d.suffix.lower() not in (".fa", ".fas", ".fasta", ".yml", ".yaml"):
                msg = (
                    f"Unable to parse filetype {d.suffix}, "
                    "please provide a .fasta or .yaml file."
                )
            else:
                valid.append(d)
                continue

            if not skip_bad_inputs:
                raise RuntimeError(msg)
            click.echo(f"Skipping: {msg}")

        data = valid
    else:
        data = [data]

    return data


def filter_inputs_structure(
    manifest: Manifest,
    outdir: Path,
    override: bool = False,
) -> Manifest:
    """Filter the manifest to only include missing predictions.

    Parameters
    ----------
    manifest : Manifest
        The manifest of the input data.
    outdir : Path
        The output directory.
    override: bool
        Whether to override existing predictions.

    Returns
    -------
    Manifest
        The manifest of the filtered input data.

    """
    from boltz.data.types import Manifest

    # Check if existing predictions are found (only top-level prediction folders)
    pred_dir = outdir / "predictions"
    if pred_dir.exists():
        existing = {d.name for d in pred_dir.iterdir() if d.is_dir()}
    else:
        existing = set()

    # Remove them from the input data, but keep records that need
    # a pre_affinity npz that wasn't written (e.g. prior run failed)
    if existing and not override:
        keep = []
        for r in manifest.records:
            if r.id not in existing:
                keep.append(r)
            elif r.affinity and not (
                pred_dir / r.id / f"pre_affinity_{r.id}.npz"
            ).exists():
                # Prior run created the directory but failed to write
                # the pre_affinity file needed for the affinity pass
                keep.append(r)
        manifest = Manifest(keep)
        msg = (
            f"Found some existing predictions ({len(existing)}), "
            f"skipping and running only the missing ones, "
            "if any. If you wish to override these existing "
            "predictions, please set the --override flag."
        )
        click.echo(msg)
    elif existing and override:
        msg = f"Found {len(existing)} existing predictions, will override."
        click.echo(msg)

    return manifest


def filter_inputs_affinity(
    manifest: Manifest,
    outdir: Path,
    override: bool = False,
) -> Manifest:
    """Check the input data and output directory for affinity.

    Parameters
    ----------
    manifest : Manifest
        The manifest.
    outdir : Path
        The output directory.
    override: bool
        Whether to override existing predictions.

    Returns
    -------
    Manifest
        The manifest of the filtered input data.

    """
    from boltz.data.types import Manifest, get_affinity_output_id

    click.echo("Checking input data for affinity.")

    pred_dir = outdir / "predictions"

    def affinity_output_path(record):
        output_id = get_affinity_output_id(record)
        return pred_dir / record.id / f"affinity_{output_id}.json"

    existing = [
        record
        for record in manifest.records
        if record.affinity and affinity_output_path(record).exists()
    ]

    # Remove them from the input data
    if existing and not override:
        existing_keys = {
            (record.id, get_affinity_output_id(record)) for record in existing
        }
        manifest = Manifest(
            [
                record
                for record in manifest.records
                if (record.id, get_affinity_output_id(record)) not in existing_keys
            ]
        )
        num_skipped = len(existing)
        msg = (
            f"Found some existing affinity predictions ({num_skipped}), "
            f"skipping and running only the missing ones, "
            "if any. If you wish to override these existing "
            "affinity predictions, please set the --override flag."
        )
        click.echo(msg)
    elif existing and override:
        msg = "Found existing affinity predictions, will override."
        click.echo(msg)

    # Skip records whose structure prediction failed (no pre_affinity npz)
    missing = [
        r for r in manifest.records
        if r.affinity
        and not (pred_dir / r.id / f"pre_affinity_{r.id}.npz").exists()
    ]
    if missing:
        missing_ids = {r.id for r in missing}
        manifest = Manifest([r for r in manifest.records if r.id not in missing_ids])
        click.echo(
            f"Warning: {len(missing)} record(s) requested affinity but have no "
            f"structure prediction (pre_affinity_*.npz missing). Skipping: "
            f"{', '.join(r.id for r in missing)}"
        )

    return manifest


def expand_affinity_records(manifest: Manifest) -> Manifest:
    """Expand repeated-binder affinity requests into per-copy records."""
    from boltz.data import const
    from boltz.data.types import AffinityInfo, Manifest

    expanded_records = []
    for record in manifest.records:
        affinity = record.affinity
        if affinity is None:
            expanded_records.append(record)
            continue

        if affinity.chain_name is not None:
            binder_chain = next(
                (
                    chain
                    for chain in record.chains
                    if chain.chain_name == affinity.chain_name
                ),
                None,
            )
            if binder_chain is None:
                msg = (
                    f"Affinity target mismatch for record {record.id}: "
                    f"could not find binder chain {affinity.chain_name!r} in "
                    "the parsed structure."
                )
                raise click.ClickException(msg)
        else:
            binder_chain = next(
                (
                    chain
                    for chain in record.chains
                    if chain.chain_id == affinity.chain_id
                ),
                None,
            )
            if binder_chain is None:
                msg = (
                    f"Affinity target mismatch for record {record.id}: "
                    f"could not find binder chain_id {affinity.chain_id} in "
                    "the parsed structure."
                )
                raise click.ClickException(msg)

        if binder_chain.entity_id is None:
            msg = (
                f"Affinity target mismatch for record {record.id}: "
                f"binder {binder_chain.chain_name!r} is missing an entity_id."
            )
            raise click.ClickException(msg)

        binder_copies = [
            chain
            for chain in record.chains
            if chain.entity_id == binder_chain.entity_id
            and chain.mol_type == const.chain_type_ids["NONPOLYMER"]
        ]
        if len(binder_copies) <= 1:
            expanded_records.append(record)
            continue

        for chain in binder_copies:
            sanitized_chain_name = re.sub(r"[^\w.-]", "_", chain.chain_name)
            expanded_records.append(
                replace(
                    record,
                    affinity=AffinityInfo(
                        chain_id=chain.chain_id,
                        mw=affinity.mw,
                        chain_name=chain.chain_name,
                        output_id=f"{record.id}_{sanitized_chain_name}",
                    ),
                )
            )

    return Manifest(expanded_records)


def compute_msa(
    data: dict[str, str],
    target_id: str,
    msa_dir: Path,
    msa_server_url: str,
    msa_pairing_strategy: str,
    msa_server_username: Optional[str] = None,
    msa_server_password: Optional[str] = None,
    api_key_header: Optional[str] = None,
    api_key_value: Optional[str] = None,
) -> None:
    """Compute the MSA for the input data.

    Parameters
    ----------
    data : dict[str, str]
        The input protein sequences.
    target_id : str
        The target id.
    msa_dir : Path
        The msa directory.
    msa_server_url : str
        The MSA server URL.
    msa_pairing_strategy : str
        The MSA pairing strategy.
    msa_server_username : str, optional
        Username for basic authentication with MSA server.
    msa_server_password : str, optional
        Password for basic authentication with MSA server.
    api_key_header : str, optional
        Custom header key for API key authentication (default: X-API-Key).
    api_key_value : str, optional
        Custom header value for API key authentication (overrides --api_key if set).

    """
    from boltz.data import const
    from boltz.data.msa.mmseqs2 import run_mmseqs2

    click.echo(f"Calling MSA server for target {target_id} with {len(data)} sequences")
    click.echo(f"MSA server URL: {msa_server_url}")
    click.echo(f"MSA pairing strategy: {msa_pairing_strategy}")

    # Construct auth headers if API key header/value is provided
    auth_headers = None
    if api_key_value:
        key = api_key_header if api_key_header else "X-API-Key"
        value = api_key_value
        auth_headers = {
            key: value
        }
        click.echo(f"Using API key authentication for MSA server (header: {key})")
    elif msa_server_username and msa_server_password:
        click.echo("Using basic authentication for MSA server")
    else:
        click.echo("No authentication provided for MSA server")
    
    if len(data) > 1:
        paired_msas = run_mmseqs2(
            list(data.values()),
            msa_dir / f"{target_id}_paired_tmp",
            use_env=True,
            use_pairing=True,
            host_url=msa_server_url,
            pairing_strategy=msa_pairing_strategy,
            msa_server_username=msa_server_username,
            msa_server_password=msa_server_password,
            auth_headers=auth_headers,
        )
    else:
        paired_msas = [""] * len(data)

    unpaired_msa = run_mmseqs2(
        list(data.values()),
        msa_dir / f"{target_id}_unpaired_tmp",
        use_env=True,
        use_pairing=False,
        host_url=msa_server_url,
        pairing_strategy=msa_pairing_strategy,
        msa_server_username=msa_server_username,
        msa_server_password=msa_server_password,
        auth_headers=auth_headers,
    )

    for idx, name in enumerate(data):
        # Get paired sequences
        paired = paired_msas[idx].strip().splitlines()
        paired = paired[1::2]  # ignore headers
        paired = paired[: const.max_paired_seqs]

        # Set key per row and remove empty sequences
        keys = [idx for idx, s in enumerate(paired) if s != "-" * len(s)]
        paired = [s for s in paired if s != "-" * len(s)]

        # Combine paired-unpaired sequences
        unpaired = unpaired_msa[idx].strip().splitlines()
        unpaired = unpaired[1::2]
        unpaired = unpaired[: (const.max_msa_seqs - len(paired))]
        if paired:
            unpaired = unpaired[1:]  # ignore query is already present

        # Combine
        seqs = paired + unpaired
        keys = keys + [-1] * len(unpaired)

        # Dump MSA
        csv_str = ["key,sequence"] + [f"{key},{seq}" for key, seq in zip(keys, seqs)]

        msa_path = msa_dir / f"{name}.csv"
        with msa_path.open("w") as f:
            f.write("\n".join(csv_str))


def process_input(  # noqa: C901, PLR0912, PLR0915, D103
    path: Path,
    ccd: dict,
    msa_dir: Path,
    mol_dir: Path,
    boltz2: bool,
    use_msa_server: bool,
    msa_server_url: str,
    msa_pairing_strategy: str,
    msa_server_username: Optional[str],
    msa_server_password: Optional[str],
    api_key_header: Optional[str],
    api_key_value: Optional[str],
    max_msa_seqs: int,
    processed_msa_dir: Path,
    processed_constraints_dir: Path,
    processed_templates_dir: Path,
    processed_mols_dir: Path,
    structure_dir: Path,
    records_dir: Path,
    skip_bad_inputs: bool = False,
) -> Optional[str]:
    try:
        # Parse data
        if path.suffix.lower() in (".fa", ".fas", ".fasta"):
            from boltz.data.parse.fasta import parse_fasta

            target = parse_fasta(path, ccd, mol_dir, boltz2)
        elif path.suffix.lower() in (".yml", ".yaml"):
            from boltz.data.parse.yaml import parse_yaml

            target = parse_yaml(path, ccd, mol_dir, boltz2)
        elif path.is_dir():
            msg = f"Found directory {path} instead of .fasta or .yaml, skipping."
            raise RuntimeError(msg)  # noqa: TRY301
        else:
            msg = (
                f"Unable to parse filetype {path.suffix}, "
                "please provide a .fasta or .yaml file."
            )
            raise RuntimeError(msg)  # noqa: TRY301

        from rdkit import Chem

        from boltz.data import const
        from boltz.data.parse.a3m import parse_a3m
        from boltz.data.parse.csv import parse_csv
        from boltz.data.types import MSA

        # Get target id
        target_id = target.record.id

        # Get all MSA ids and decide whether to generate MSA
        to_generate = {}
        prot_id = const.chain_type_ids["PROTEIN"]
        for chain in target.record.chains:
            # Add to generate list, assigning entity id
            if (chain.mol_type == prot_id) and (chain.msa_id == 0):
                entity_id = chain.entity_id
                msa_id = f"{target_id}_{entity_id}"
                to_generate[msa_id] = target.sequences[entity_id]
                chain.msa_id = msa_dir / f"{msa_id}.csv"

            # We do not support msa generation for non-protein chains
            elif chain.msa_id == 0:
                chain.msa_id = -1

        # Generate MSA
        if to_generate and not use_msa_server:
            msg = "Missing MSA's in input and --use_msa_server flag not set."
            raise RuntimeError(msg)  # noqa: TRY301

        if to_generate:
            msg = f"Generating MSA for {path} with {len(to_generate)} protein entities."
            click.echo(msg)
            compute_msa(
                data=to_generate,
                target_id=target_id,
                msa_dir=msa_dir,
                msa_server_url=msa_server_url,
                msa_pairing_strategy=msa_pairing_strategy,
                msa_server_username=msa_server_username,
                msa_server_password=msa_server_password,
                api_key_header=api_key_header,
                api_key_value=api_key_value,
            )

        # Parse MSA data
        msas = sorted({c.msa_id for c in target.record.chains if c.msa_id != -1})
        msa_id_map = {}
        for msa_idx, msa_id in enumerate(msas):
            # Check that raw MSA exists
            msa_path = Path(msa_id)
            if not msa_path.exists():
                msg = f"MSA file {msa_path} not found."
                raise FileNotFoundError(msg)  # noqa: TRY301

            # Dump processed MSA
            processed = processed_msa_dir / f"{target_id}_{msa_idx}.npz"
            msa_id_map[msa_id] = f"{target_id}_{msa_idx}"
            if not processed.exists():
                # Parse A3M
                if msa_path.suffix == ".a3m":
                    msa: MSA = parse_a3m(
                        msa_path,
                        taxonomy=None,
                        max_seqs=max_msa_seqs,
                    )
                elif msa_path.suffix == ".csv":
                    msa: MSA = parse_csv(msa_path, max_seqs=max_msa_seqs)
                else:
                    msg = f"MSA file {msa_path} not supported, only a3m or csv."
                    raise RuntimeError(msg)  # noqa: TRY301

                msa.dump(processed)

        # Modify records to point to processed MSA
        for c in target.record.chains:
            if (c.msa_id != -1) and (c.msa_id in msa_id_map):
                c.msa_id = msa_id_map[c.msa_id]

        # Dump templates
        for template_id, template in target.templates.items():
            name = f"{target.record.id}_{template_id}.npz"
            template_path = processed_templates_dir / name
            template.dump(template_path)

        # Dump constraints
        constraints_path = processed_constraints_dir / f"{target.record.id}.npz"
        target.residue_constraints.dump(constraints_path)

        # Dump extra molecules
        Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
        with (processed_mols_dir / f"{target.record.id}.pkl").open("wb") as f:
            pickle.dump(target.extra_mols, f)

        # Dump structure
        struct_path = structure_dir / f"{target.record.id}.npz"
        target.structure.dump(struct_path)

        # Dump record
        record_path = records_dir / f"{target.record.id}.json"
        target.record.dump(record_path)

        return None

    except Exception as e:  # noqa: BLE001
        if not skip_bad_inputs:
            raise
        return f"{path.name}: {e}"


@rank_zero_only
def process_inputs(
    data: list[Path],
    out_dir: Path,
    ccd_path: Path,
    mol_dir: Path,
    msa_server_url: str,
    msa_pairing_strategy: str,
    max_msa_seqs: int = 8192,
    use_msa_server: bool = False,
    msa_server_username: Optional[str] = None,
    msa_server_password: Optional[str] = None,
    api_key_header: Optional[str] = None,
    api_key_value: Optional[str] = None,
    boltz2: bool = False,
    preprocessing_threads: int = 1,
    skip_bad_inputs: bool = False,
) -> Manifest:
    """Process the input data and output directory.

    Parameters
    ----------
    data : list[Path]
        The input data.
    out_dir : Path
        The output directory.
    ccd_path : Path
        The path to the CCD dictionary.
    max_msa_seqs : int, optional
        Max number of MSA sequences, by default 8192.
    use_msa_server : bool, optional
        Whether to use the MMSeqs2 server for MSA generation, by default False.
    msa_server_username : str, optional
        Username for basic authentication with MSA server, by default None.
    msa_server_password : str, optional
        Password for basic authentication with MSA server, by default None.
    api_key_header : str, optional
        Custom header key for API key authentication (default: X-API-Key).
    api_key_value : str, optional
        Custom header value for API key authentication (overrides --api_key if set).
    boltz2: bool, optional
        Whether to use Boltz2, by default False.
    preprocessing_threads: int, optional
        The number of threads to use for preprocessing, by default 1.

    Returns
    -------
    Manifest
        The manifest of the processed input data.

    """
    from tqdm import tqdm

    from boltz.data.mol import load_canonicals
    from boltz.data.types import Manifest, Record

    # Validate mutually exclusive authentication methods
    has_basic_auth = msa_server_username and msa_server_password
    has_api_key = api_key_value is not None
    
    if has_basic_auth and has_api_key:
        raise ValueError(
            "Cannot use both basic authentication (--msa_server_username/--msa_server_password) "
            "and API key authentication (--api_key_header/--api_key_value). Please use only one authentication method."
        )

    # Check if records exist at output path
    records_dir = out_dir / "processed" / "records"
    if records_dir.exists():
        # Load existing records
        existing = [Record.load(p) for p in records_dir.glob("*.json")]
        processed_ids = {record.id for record in existing}

        # Filter to missing only
        data = [d for d in data if d.stem not in processed_ids]

        # Nothing to do, update the manifest and return
        if data:
            click.echo(
                f"Found {len(existing)} existing processed inputs, skipping them."
            )
        else:
            click.echo("All inputs are already processed.")
            updated_manifest = Manifest(existing)
            updated_manifest.dump(out_dir / "processed" / "manifest.json")

    # Create output directories
    msa_dir = out_dir / "msa"
    records_dir = out_dir / "processed" / "records"
    structure_dir = out_dir / "processed" / "structures"
    processed_msa_dir = out_dir / "processed" / "msa"
    processed_constraints_dir = out_dir / "processed" / "constraints"
    processed_templates_dir = out_dir / "processed" / "templates"
    processed_mols_dir = out_dir / "processed" / "mols"
    predictions_dir = out_dir / "predictions"

    out_dir.mkdir(parents=True, exist_ok=True)
    msa_dir.mkdir(parents=True, exist_ok=True)
    records_dir.mkdir(parents=True, exist_ok=True)
    structure_dir.mkdir(parents=True, exist_ok=True)
    processed_msa_dir.mkdir(parents=True, exist_ok=True)
    processed_constraints_dir.mkdir(parents=True, exist_ok=True)
    processed_templates_dir.mkdir(parents=True, exist_ok=True)
    processed_mols_dir.mkdir(parents=True, exist_ok=True)
    predictions_dir.mkdir(parents=True, exist_ok=True)

    # Load CCD
    if boltz2:
        ccd = load_canonicals(mol_dir)
    else:
        with ccd_path.open("rb") as file:
            ccd = pickle.load(file)  # noqa: S301

    # Create partial function
    process_input_partial = partial(
        process_input,
        ccd=ccd,
        msa_dir=msa_dir,
        mol_dir=mol_dir,
        boltz2=boltz2,
        use_msa_server=use_msa_server,
        msa_server_url=msa_server_url,
        msa_pairing_strategy=msa_pairing_strategy,
        msa_server_username=msa_server_username,
        msa_server_password=msa_server_password,
        api_key_header=api_key_header,
        api_key_value=api_key_value,
        max_msa_seqs=max_msa_seqs,
        processed_msa_dir=processed_msa_dir,
        processed_constraints_dir=processed_constraints_dir,
        processed_templates_dir=processed_templates_dir,
        processed_mols_dir=processed_mols_dir,
        structure_dir=structure_dir,
        records_dir=records_dir,
        skip_bad_inputs=skip_bad_inputs,
    )

    # Parse input data
    preprocessing_threads = min(preprocessing_threads, len(data))
    click.echo(f"Processing {len(data)} inputs with {preprocessing_threads} threads.")

    if preprocessing_threads > 1 and len(data) > 1:
        with Pool(preprocessing_threads) as pool:
            errors = list(tqdm(pool.imap(process_input_partial, data), total=len(data)))
    else:
        errors = []
        for path in tqdm(data):
            errors.append(process_input_partial(path))

    # Report failures
    failures = [e for e in errors if e is not None]
    if failures:
        click.echo(f"\n{len(failures)} input(s) failed:")
        for msg in failures:
            click.echo(f"  {msg}")

    # Load all records and write manifest
    records = [Record.load(p) for p in records_dir.glob("*.json")]
    manifest = Manifest(records)
    manifest.dump(out_dir / "processed" / "manifest.json")


@click.group()
@click.version_option(version=_pkg_version("boltz-community"), prog_name="boltz-community")
def cli() -> None:
    """Boltz."""
    return


def _parse_devices(value: str) -> Union[int, List[int]]:
    """Parse --devices as an integer count or comma-separated list of device IDs."""
    if "," in value:
        return [int(x.strip()) for x in value.split(",")]
    return int(value)


@cli.command(context_settings={"show_default": True})
@click.argument("data", type=click.Path(exists=True))
@click.option(
    "--out_dir",
    type=click.Path(exists=False),
    help="The path where to save the predictions.",
    default="./",
)
@click.option(
    "--cache",
    type=click.Path(exists=False),
    help=(
        "The directory where to download the data and model. "
        "Default is ~/.boltz, or $BOLTZ_CACHE if set."
    ),
    default=get_cache_path,
)
@click.option(
    "--checkpoint",
    type=click.Path(exists=True),
    help="An optional checkpoint, will use the provided Boltz-1 model by default.",
    default=None,
)
@click.option(
    "--devices",
    type=str,
    callback=lambda ctx, param, value: _parse_devices(value),
    help=(
        "Devices to use for prediction. Pass an integer to set the number of devices "
        "(e.g. --devices 2 uses 2 GPUs), or a comma-separated list of specific device IDs "
        "(e.g. --devices 0,1 uses GPUs 0 and 1). To target a single GPU by index, use "
        "CUDA_VISIBLE_DEVICES=1 boltz predict ... Default is 1."
    ),
    default="1",
)
@click.option(
    "--accelerator",
    type=click.Choice(["gpu", "cpu", "tpu", "mps"]),
    help="The accelerator to use for prediction.",
    default="gpu",
)
@click.option(
    "--recycling_steps",
    type=int,
    help="The number of recycling steps to use for prediction.",
    default=3,
)
@click.option(
    "--sampling_steps",
    type=int,
    help="The number of sampling steps to use for prediction.",
    default=200,
)
@click.option(
    "--diffusion_samples",
    type=int,
    help="The number of diffusion samples to use for prediction.",
    default=1,
)
@click.option(
    "--max_parallel_samples",
    type=int,
    help="The maximum number of samples to predict in parallel.",
    default=5,
)
@click.option(
    "--step_scale",
    type=float,
    help=(
        "The step size is related to the temperature at "
        "which the diffusion process samples the distribution. "
        "The lower the higher the diversity among samples "
        "(recommended between 1 and 2). "
        "Default is 1.638 for Boltz-1 and 1.5 for Boltz-2. "
        "If not provided, the default step size will be used."
    ),
    default=None,
)
@click.option(
    "--write_full_pae/--no_write_full_pae",
    type=bool,
    is_flag=True,
    default=True,
    help="Whether to dump the pae into a npz file.",
)
@click.option(
    "--write_full_pde/--no_write_full_pde",
    type=bool,
    is_flag=True,
    default=True,
    help="Whether to dump the pde into a npz file.",
)
@click.option(
    "--output_format",
    type=click.Choice(["pdb", "mmcif"]),
    help="The output format to use for the predictions.",
    default="mmcif",
)
@click.option(
    "--num_workers",
    type=int,
    help="The number of dataloader workers to use for prediction.",
    default=2,
)
@click.option(
    "--batch_size",
    type=click.IntRange(min=1),
    help=(
        "The number of inputs to process per batch during Boltz-2 structure "
        "prediction. Affinity prediction remains batch_size=1."
    ),
    default=1,
)
@click.option(
    "--override",
    is_flag=True,
    help="Override existing found predictions.",
)
@click.option(
    "--skip_bad_inputs",
    is_flag=True,
    help="Skip invalid inputs instead of aborting. Default: abort on first error.",
)
@click.option(
    "--seed",
    type=int,
    help="Seed to use for random number generator.",
    default=None,
)
@click.option(
    "--use_msa_server",
    is_flag=True,
    help="Whether to use the MMSeqs2 server for MSA generation.",
)
@click.option(
    "--msa_server_url",
    type=str,
    help="MSA server url. Used only if --use_msa_server is set. ",
    default="https://api.colabfold.com",
)
@click.option(
    "--msa_pairing_strategy",
    type=click.Choice(["greedy", "complete"]),
    help="Pairing strategy to use. Used only if --use_msa_server is set.",
    default="greedy",
)
@click.option(
    "--msa_server_username",
    type=str,
    help="MSA server username for basic auth. Used only if --use_msa_server is set. Can also be set via BOLTZ_MSA_USERNAME environment variable.",
    default=None,
)
@click.option(
    "--msa_server_password",
    type=str,
    help="MSA server password for basic auth. Used only if --use_msa_server is set. Can also be set via BOLTZ_MSA_PASSWORD environment variable.",
    default=None,
)
@click.option(
    "--api_key_header",
    type=str,
    help="Custom header key for API key authentication (default: X-API-Key).",
    default=None,
)
@click.option(
    "--api_key_value",
    type=str,
    help="Custom header value for API key authentication.",
    default=None,
)
@click.option(
    "--use_potentials",
    is_flag=True,
    help="Whether to use potentials for steering.",
)
@click.option(
    "--model",
    default="boltz2",
    type=click.Choice(["boltz1", "boltz2"]),
    help="The model to use for prediction.",
)
@click.option(
    "--method",
    type=str,
    help="The method to use for prediction.",
    default=None,
)
@click.option(
    "--preprocessing-threads",
    type=int,
    help="The number of threads to use for preprocessing.",
    default=_available_cpu_count(),
)
@click.option(
    "--affinity_mw_correction",
    is_flag=True,
    type=bool,
    help="Whether to add the Molecular Weight correction to the affinity value head.",
)
@click.option(
    "--sampling_steps_affinity",
    type=int,
    help="The number of sampling steps to use for affinity prediction.",
    default=200,
)
@click.option(
    "--diffusion_samples_affinity",
    type=int,
    help="The number of diffusion samples to use for affinity prediction.",
    default=5,
)
@click.option(
    "--affinity_checkpoint",
    type=click.Path(exists=True),
    help="An optional checkpoint, will use the provided Boltz-1 model by default.",
    default=None,
)
@click.option(
    "--max_msa_seqs",
    type=int,
    help="The maximum number of MSA sequences to use for prediction.",
    default=8192,
)
@click.option(
    "--subsample_msa/--no_subsample_msa",
    default=True,
    help="Whether to subsample the MSA.",
)
@click.option(
    "--num_subsampled_msa",
    type=int,
    help="The number of MSA sequences to subsample.",
    default=1024,
)
@click.option(
    "--no_kernels",
    is_flag=True,
    help="Whether to disable the kernels.",
)
@click.option(
    "--flash_attn",
    is_flag=True,
    help=(
        "Enable PyTorch scaled_dot_product_attention (FlashAttention-2 / "
        "memory-efficient attention) for all attention layers. "
        "Requires an Ampere (sm_80+) or newer GPU and PyTorch ≥ 2.0. "
        "Provides O(N) memory usage vs O(N²) and significant throughput gains "
        "on large sequences. Default False."
    ),
)
@click.option(
    "--write_embeddings",
    is_flag=True,
    help="Whether to dump the s and z embeddings into a npz file.",
)
def predict(  # noqa: C901, PLR0915, PLR0912
    data: str,
    out_dir: str,
    cache: str,
    checkpoint: Optional[str],
    affinity_checkpoint: Optional[str],
    devices: Union[int, List[int]],
    accelerator: str,
    recycling_steps: int,
    sampling_steps: int,
    diffusion_samples: int,
    sampling_steps_affinity: int,
    diffusion_samples_affinity: int,
    max_parallel_samples: Optional[int],
    step_scale: Optional[float],
    write_full_pae: bool,
    write_full_pde: bool,
    output_format: Literal["pdb", "mmcif"],
    num_workers: int,
    batch_size: int,
    override: bool,
    skip_bad_inputs: bool,
    seed: Optional[int],
    use_msa_server: bool,
    msa_server_url: str,
    msa_pairing_strategy: str,
    msa_server_username: Optional[str],
    msa_server_password: Optional[str],
    api_key_header: Optional[str],
    api_key_value: Optional[str],
    use_potentials: bool,
    model: Literal["boltz1", "boltz2"],
    method: Optional[str],
    affinity_mw_correction: Optional[bool],
    preprocessing_threads: int,
    max_msa_seqs: int,
    subsample_msa: bool,
    num_subsampled_msa: int,
    no_kernels: bool,
    flash_attn: bool,
    write_embeddings: bool,
) -> None:
    """Run predictions with Boltz."""
    import torch
    from pytorch_lightning import Trainer, seed_everything
    from pytorch_lightning.strategies import DDPStrategy
    from rdkit import Chem

    from boltz.data import const
    from boltz.data.module.inference import BoltzInferenceDataModule
    from boltz.data.module.inferencev2 import Boltz2InferenceDataModule
    from boltz.data.types import Manifest
    from boltz.data.write.writer import BoltzAffinityWriter, BoltzWriter
    from boltz.model.models.boltz1 import Boltz1
    from boltz.model.models.boltz2 import Boltz2

    accelerator = _normalize_accelerator(accelerator)

    # PyTorch 2.6+ defaults torch.load to weights_only=True, which rejects
    # Lightning checkpoints containing OmegaConf config objects and custom types.
    # Since boltz only loads its own checkpoints from trusted sources, restore
    # the pre-2.6 default.
    if not getattr(torch.load, "_boltz_patched", False):
        _original_torch_load = torch.load

        def _patched_load(*a, _orig=_original_torch_load, **kw):
            location = kw.get("map_location")
            try:
                location_text = str(location)
            except Exception:
                location_text = ""

            if torch.backends.mps.is_available() and not hasattr(torch.mps, "current_device"):
                if location_text.startswith(("mps", "cuda")):
                    kw["map_location"] = "cpu"
            return _orig(*a, **{**kw, "weights_only": False})

        _patched_load._boltz_patched = True  # type: ignore[attr-defined]
        torch.load = _patched_load  # type: ignore[assignment]

    # Some dev/nightly MPS builds expose Metal support but omit
    # torch.mps.current_device(), which breaks checkpoint deserialization.
    # Apply the guard unconditionally when MPS is available so both explicit
    # --accelerator mps and default GPU runs stay compatible.
    if torch.backends.mps.is_available():
        _ensure_torch_mps_current_device_compat()

    # If cpu, write a friendly warning
    if accelerator == "cpu":
        msg = "Running on CPU, this will be slow. Consider using a GPU."
        click.echo(msg)

    if batch_size > 1 and model != "boltz2":
        raise click.UsageError(
            "Batched structure inference is only supported for Boltz-2."
        )
    if batch_size > 1 and use_potentials:
        raise click.UsageError(
            "Batched structure inference with --use_potentials is not supported yet."
        )
    if batch_size > 1:
        click.echo(
            "Warning: contact guidance is not yet supported with batched "
            "Boltz-2 structure inference; running without it. Use "
            "--batch_size 1 for guided inference."
        )

    # MPS (Apple Silicon) validation and info
    if accelerator == "mps":
        if not torch.backends.mps.is_available():
            msg = (
                "MPS accelerator requested but not available. "
                "MPS requires macOS 12.3+ on Apple Silicon."
            )
            raise click.UsageError(msg)
        click.echo("Running on Apple Silicon GPU (MPS).")

    # Supress some lightning warnings
    warnings.filterwarnings(
        "ignore", ".*that has Tensor Cores. To properly utilize them.*"
    )

    # Set no grad
    torch.set_grad_enabled(False)

    # Ignore matmul precision warning
    torch.set_float32_matmul_precision("highest")

    # Set rdkit pickle logic
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    # Set seed if desired
    if seed is not None:
        seed_everything(seed)

    for key in ["CUEQ_DEFAULT_CONFIG", "CUEQ_DISABLE_AOT_TUNING"]:
        # Disable kernel tuning by default,
        # but do not modify envvar if already set by caller
        os.environ[key] = os.environ.get(key, "1")

    # Set cache path
    cache = Path(cache).expanduser()
    cache.mkdir(parents=True, exist_ok=True)

    # Get MSA server credentials from environment variables if not provided
    if use_msa_server:
        if msa_server_username is None:
            msa_server_username = os.environ.get("BOLTZ_MSA_USERNAME")
        if msa_server_password is None:
            msa_server_password = os.environ.get("BOLTZ_MSA_PASSWORD")
        if api_key_value is None:
            api_key_value = os.environ.get("MSA_API_KEY_VALUE")
        
        click.echo(f"MSA server enabled: {msa_server_url}")
        if api_key_value:
            click.echo("MSA server authentication: using API key header")
        elif msa_server_username and msa_server_password:
            click.echo("MSA server authentication: using basic auth")
        else:
            click.echo("MSA server authentication: no credentials provided")

    # Create output directories
    data = Path(data).expanduser()
    out_dir = Path(out_dir).expanduser()
    out_dir = out_dir / f"boltz_results_{data.stem}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Download necessary data and model
    if model == "boltz1":
        download_boltz1(cache)
    elif model == "boltz2":
        download_boltz2(cache)
    else:
        msg = f"Model {model} not supported. Supported: boltz1, boltz2."
        raise ValueError(f"Model {model} not supported.")

    # Validate inputs
    data = check_inputs(data, skip_bad_inputs=skip_bad_inputs)

    # Check method
    if method is not None:
        if model == "boltz1":
            msg = "Method conditioning is not supported for Boltz-1."
            raise ValueError(msg)
        if method.lower() not in const.method_types_ids:
            method_names = list(const.method_types_ids.keys())
            msg = f"Method {method} not supported. Supported: {method_names}"
            raise ValueError(msg)

    # Process inputs
    ccd_path = cache / "ccd.pkl"
    mol_dir = cache / "mols"
    process_inputs(
        data=data,
        out_dir=out_dir,
        ccd_path=ccd_path,
        mol_dir=mol_dir,
        use_msa_server=use_msa_server,
        msa_server_url=msa_server_url,
        msa_pairing_strategy=msa_pairing_strategy,
        msa_server_username=msa_server_username,
        msa_server_password=msa_server_password,
        api_key_header=api_key_header,
        api_key_value=api_key_value,
        boltz2=model == "boltz2",
        preprocessing_threads=preprocessing_threads,
        max_msa_seqs=max_msa_seqs,
        skip_bad_inputs=skip_bad_inputs,
    )

    # Load manifest
    manifest = Manifest.load(out_dir / "processed" / "manifest.json")

    if not manifest.records:
        raise click.ClickException(
            "No inputs were successfully processed. "
            "Check the error messages above for details."
        )

    # Filter out existing predictions
    filtered_manifest = filter_inputs_structure(
        manifest=manifest,
        outdir=out_dir,
        override=override,
    )

    # Load processed data
    processed_dir = out_dir / "processed"
    processed = BoltzProcessedInput(
        manifest=filtered_manifest,
        targets_dir=processed_dir / "structures",
        msa_dir=processed_dir / "msa",
        constraints_dir=(
            (processed_dir / "constraints")
            if (processed_dir / "constraints").exists()
            else None
        ),
        template_dir=(
            (processed_dir / "templates")
            if (processed_dir / "templates").exists()
            else None
        ),
        extra_mols_dir=(
            (processed_dir / "mols") if (processed_dir / "mols").exists() else None
        ),
    )

    # Set up trainer
    strategy = "auto"
    if accelerator == "mps":
        # MPS only supports single-device, float32
        if (isinstance(devices, int) and devices > 1) or (
            isinstance(devices, list) and len(devices) > 1
        ):
            click.echo("Warning: MPS only supports a single device, ignoring --devices.")
        devices = 1
        strategy = "auto"
    elif (isinstance(devices, int) and devices > 1) or (
        isinstance(devices, list) and len(devices) > 1
    ):
        start_method = "fork" if platform.system() != "win32" and platform.system() != "Windows" else "spawn"
        strategy = DDPStrategy(start_method=start_method)
        device_count = len(devices) if isinstance(devices, list) else devices
        if len(filtered_manifest.records) < device_count:
            msg = (
                "Number of requested devices is greater "
                "than the number of predictions, taking the minimum."
            )
            click.echo(msg)
            if isinstance(devices, list):
                devices = devices[: max(1, len(filtered_manifest.records))]
            else:
                devices = max(1, min(len(filtered_manifest.records), devices))

    # Set up model parameters
    if model == "boltz2":
        diffusion_params = Boltz2DiffusionParams()
        step_scale = 1.5 if step_scale is None else step_scale
        diffusion_params.step_scale = step_scale
        pairformer_args = PairformerArgsV2()
    else:
        diffusion_params = BoltzDiffusionParams()
        step_scale = 1.638 if step_scale is None else step_scale
        diffusion_params.step_scale = step_scale
        pairformer_args = PairformerArgs()

    msa_args = MSAModuleArgs(
        subsample_msa=subsample_msa,
        num_subsampled_msa=num_subsampled_msa,
        use_paired_feature=model == "boltz2",
    )

    # Create prediction writer
    pred_writer = BoltzWriter(
        data_dir=processed.targets_dir,
        output_dir=out_dir / "predictions",
        output_format=output_format,
        boltz2=model == "boltz2",
        write_embeddings=write_embeddings,
    )

    # MPS and CPU do not support bf16-mixed
    if accelerator in ("mps", "cpu"):
        precision = 32
    elif model == "boltz1":
        precision = 32
    else:
        precision = "bf16-mixed"

    # Set up trainer
    trainer = Trainer(
        default_root_dir=out_dir,
        strategy=strategy,
        callbacks=[pred_writer],
        accelerator=accelerator,
        devices=devices,
        precision=precision,
    )

    map_location = _resolve_checkpoint_load_map_location(accelerator)

    pin_memory = accelerator not in ("cpu", "mps")

    if filtered_manifest.records:
        msg = f"Running structure prediction for {len(filtered_manifest.records)} input"
        msg += "s." if len(filtered_manifest.records) > 1 else "."
        click.echo(msg)

        # Create data module
        if model == "boltz2":
            data_module = Boltz2InferenceDataModule(
                manifest=processed.manifest,
                target_dir=processed.targets_dir,
                msa_dir=processed.msa_dir,
                mol_dir=mol_dir,
                num_workers=num_workers,
                batch_size=batch_size,
                constraints_dir=processed.constraints_dir,
                template_dir=processed.template_dir,
                extra_mols_dir=processed.extra_mols_dir,
                override_method=method,
                pin_memory=pin_memory,
            )
        else:
            data_module = BoltzInferenceDataModule(
                manifest=processed.manifest,
                target_dir=processed.targets_dir,
                msa_dir=processed.msa_dir,
                num_workers=num_workers,
                constraints_dir=processed.constraints_dir,
                pin_memory=pin_memory,
            )

        # Load model
        if checkpoint is None:
            if model == "boltz2":
                checkpoint = cache / "boltz2_conf.ckpt"
            else:
                checkpoint = cache / "boltz1_conf.ckpt"
        checkpoint = Path(checkpoint).expanduser()

        predict_args = {
            "recycling_steps": recycling_steps,
            "sampling_steps": sampling_steps,
            "diffusion_samples": diffusion_samples,
            "max_parallel_samples": max_parallel_samples,
            "write_confidence_summary": True,
            "write_full_pae": write_full_pae,
            "write_full_pde": write_full_pde,
        }

        steering_args = BoltzSteeringParams()
        steering_args.fk_steering = use_potentials
        steering_args.physical_guidance_update = use_potentials
        # Contact guidance assumes per-record sequential inference today, so
        # keep it disabled for batched structure prediction until that path is audited.
        steering_args.contact_guidance_update = batch_size == 1

        model_cls = Boltz2 if model == "boltz2" else Boltz1
        if model == "boltz2":
            _validate_checkpoint_for_load(checkpoint, "Boltz-2 weights")
        else:
            _validate_checkpoint_for_load(checkpoint, "Boltz-1 weights")
        model_module = model_cls.load_from_checkpoint(
            checkpoint,
            strict=True,
            predict_args=predict_args,
            map_location=map_location,
            diffusion_process_args=asdict(diffusion_params),
            ema=False,
            use_kernels=not no_kernels,
            use_flash_attn=flash_attn,
            pairformer_args=asdict(pairformer_args),
            msa_args=asdict(msa_args),
            steering_args=asdict(steering_args),
        )
        model_module.eval()

        # Compute structure predictions
        trainer.predict(
            model_module,
            datamodule=data_module,
            return_predictions=False,
        )

        if pred_writer.failed:
            click.echo(
                f"\nWarning: {pred_writer.failed} structure prediction(s) failed "
                f"(e.g. out of memory). Affinity predictions for those records "
                f"will be skipped."
            )

    # Check if affinity predictions are needed
    if any(r.affinity for r in manifest.records):
        # Print header
        click.echo("\nPredicting property: affinity\n")

        manifest_affinity = expand_affinity_records(manifest)

        # Validate inputs
        manifest_filtered = filter_inputs_affinity(
            manifest=manifest_affinity,
            outdir=out_dir,
            override=override,
        )
        if not manifest_filtered.records:
            click.echo("Found existing affinity predictions for all inputs, skipping.")
            return

        msg = (
            f"Running affinity prediction for "
            f"{len(manifest_filtered.records)} affinity target"
        )
        msg += "s." if len(manifest_filtered.records) > 1 else "."
        click.echo(msg)

        pred_writer = BoltzAffinityWriter(
            data_dir=processed.targets_dir,
            output_dir=out_dir / "predictions",
        )

        data_module = Boltz2InferenceDataModule(
            manifest=manifest_filtered,
            target_dir=out_dir / "predictions",
            msa_dir=processed.msa_dir,
            mol_dir=mol_dir,
            num_workers=num_workers,
            # Keep affinity single-record until the batched affinity path is audited.
            batch_size=1,
            constraints_dir=processed.constraints_dir,
            template_dir=processed.template_dir,
            extra_mols_dir=processed.extra_mols_dir,
            override_method="other",
            affinity=True,
            pin_memory=pin_memory,
        )

        predict_affinity_args = {
            "recycling_steps": 5,
            "sampling_steps": sampling_steps_affinity,
            "diffusion_samples": diffusion_samples_affinity,
            "max_parallel_samples": _resolve_affinity_max_parallel_samples(
                max_parallel_samples,
                diffusion_samples_affinity,
            ),
            "write_confidence_summary": False,
            "write_full_pae": False,
            "write_full_pde": False,
        }

        # Load affinity model
        if affinity_checkpoint is None:
            affinity_checkpoint = cache / "boltz2_aff.ckpt"
        affinity_checkpoint = Path(affinity_checkpoint).expanduser()

        steering_args = BoltzSteeringParams()
        steering_args.fk_steering = False
        steering_args.physical_guidance_update = False
        steering_args.contact_guidance_update = False

        _validate_checkpoint_for_load(
            affinity_checkpoint, "Boltz-2 affinity weights"
        )
        model_module = Boltz2.load_from_checkpoint(
            affinity_checkpoint,
            strict=True,
            predict_args=predict_affinity_args,
            map_location=map_location,
            diffusion_process_args=asdict(diffusion_params),
            ema=False,
            pairformer_args=asdict(pairformer_args),
            msa_args=asdict(msa_args),
            steering_args=asdict(steering_args),
            affinity_mw_correction=affinity_mw_correction,
        )
        model_module.eval()

        trainer.callbacks[0] = pred_writer
        trainer.predict(
            model_module,
            datamodule=data_module,
            return_predictions=False,
        )


if __name__ == "__main__":
    cli()
