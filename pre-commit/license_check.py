# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


from __future__ import annotations

import ast
import logging
import subprocess
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from typing import Callable

from dataclasses import dataclass
from functools import partial
from pathlib import Path

import click

logger = logging.getLogger(__name__)

__all__: Sequence[str] = (
    # main license check functionality: per file & per directory
    #   (recursively, *.py filter)
    "license_check",
    "check_license_project_files",
    "Checked",
    # error types
    "LicenseCheckError",
    "HeaderNotFoundError",
    # functions that implement license checking behavior
    "append_license_header",
    "is_valid_python",
    "has_header",
    "ensure_license_starts_with_pound",
    "remove_existing_license_header",
    # default license header
    "LICENSE_HEADER",
    # to run main CLI program logic, w/o Click runner
    "main",
)

LICENSE_HEADER: str = (Path(__file__).parent / "license_header").read_text().strip()


@dataclass(frozen=True)
class HeaderNotFoundError(ValueError):
    """Error that indicates the pointed-to file does not have a valid license header."""

    pyfile: Path

    def __str__(self) -> str:  # noqa: D105
        return f"{self.pyfile.name} does not have the license header!"


LicenseCheckError = IOError | SyntaxError | HeaderNotFoundError
"""Errors that can be encountered during the license check process.

Specific errors and their underlying causes:
  - IOError: problem reading file
  - SyntaxError: the input file for license checking is not valid Python
  - HeaderNotFound: the input file was valid Python, but did not have the right
    license @ the header
"""


def license_check(  # noqa: PLR0911
    pyfile: Path,
    *,
    license_header: str = LICENSE_HEADER,
    modify: bool,
    replace: bool = False,
) -> LicenseCheckError | None:
    """Check file for license header, returning nothing on success or an error."""
    if not pyfile.is_file():
        return OSError(f"{pyfile.name} file does not exist!")

    with open(str(pyfile)) as rt:  # noqa: PTH123
        pyfile_contents: str = rt.read()

    maybe_err = is_valid_python(pyfile_contents)
    if maybe_err is not None:
        return maybe_err

    if has_header(pyfile_contents, license_header=license_header):
        return None

    if modify:
        # `pyfile` doesn't start with `license_header` text

        if replace:
            # does it start with some other license header?
            # if so, then we delete that before appending our new `license_header` text
            pyfile_contents = remove_existing_license_header(pyfile_contents)
            maybe_err = is_valid_python(pyfile_contents)
            if maybe_err is not None:
                return maybe_err

        pyfile_contents = append_license_header(pyfile_contents, license_header=license_header)
        maybe_err = is_valid_python(pyfile_contents)
        if maybe_err is not None:
            return maybe_err

        with open(str(pyfile), "w") as wt:  # noqa: PTH123
            wt.write(pyfile_contents)
        return None
    return HeaderNotFoundError(pyfile)


def is_valid_python(pyfile_contents: str) -> SyntaxError | None:
    """Validates python code. Returns None if it is valid, otherwise a SyntaxError."""
    try:
        _ = ast.parse(pyfile_contents)
    except SyntaxError as error:
        return error
    else:
        return None


def has_header(pyfile_contents: str, *, license_header: str = LICENSE_HEADER) -> bool:
    """True if the :param:`pyfile_contents` starts with the :param:`license_header`."""
    return pyfile_contents.startswith(license_header)


def append_license_header(pyfile_contents: str, *, license_header: str = LICENSE_HEADER, n_sep_lines: int = 2) -> str:
    """Appends the :param:`license_header` to the beginning of the input Python code.

    Inserts :param:`n_sep_lines` newlines between the license header & Python file
      content.
    """
    spacer = "\n" * n_sep_lines
    return f"{license_header}{spacer}{pyfile_contents}"


def remove_existing_license_header(pyfile_contents: str) -> str:
    """Heuristically removes the license header from a Python file's contents.

    Assumes that a license header is identified by a span of commented-out lines
    from the beginning of the file.  I.e. a big initial block of lines starting
    with "#" ==> a license header.

    Will always return the input without this "license header" block.
    """
    if not pyfile_contents.startswith("#") or len(pyfile_contents) == 0:
        return pyfile_contents
    lines: list[str] = pyfile_contents.split("\n")
    non_header_lines = lines[_last_index_of_header_comment_line(lines) + 1 :]
    return "\n".join(non_header_lines)


def _last_index_of_header_comment_line(lines: list[str]) -> int:
    """Return index into `lines` with the first line that doesn't start as a comment."""
    if len(lines) == 0:
        raise ValueError
    last_index_of_line_that_started_with_hash_from_beginning: int = -1
    for i, line in enumerate(lines):
        if line.startswith("#"):
            last_index_of_line_that_started_with_hash_from_beginning = i
        else:
            break
    if last_index_of_line_that_started_with_hash_from_beginning < 0:
        raise ValueError("Must supply non-empty lines of Python!")  # noqa: TRY003, EM101
    return last_index_of_line_that_started_with_hash_from_beginning


@dataclass(frozen=True)
class Checked:
    """Result of running license check across a collection of Python files."""

    noncompliant_files: Mapping[Path, LicenseCheckError]
    """Files that either don't have a license header for some reason or another.
    """

    n_files: int
    """Total number of Python files checked.
    """


def check_license_project_files(
    python_package_directory: Path, *, license_header: str, modify: bool, replace: bool
) -> Checked:
    """Check all Python files in a given directory tree, returning non-compliant files.

    Each returned file will be associated with the specific :class:`LicenseCheckError`.
    For more details,
    see :func:`license_check`.
    """
    assert python_package_directory.is_dir(), (
        "Input must be a directory of Python files, not a " f"directory: {python_package_directory}"
    )
    noncompliant_files = {}
    n_files = 0
    for pyfile in python_package_directory.rglob("*.py"):
        n_files += 1
        maybe_error = license_check(pyfile, license_header=license_header, modify=modify, replace=replace)
        if maybe_error is not None:
            noncompliant_files[pyfile] = maybe_error
    return Checked(noncompliant_files, n_files)


def ensure_license_starts_with_pound(license_header_contents: str) -> str:
    """Ensures lines of the license headers start with "# "; add if necessary."""
    if len(license_header_contents) == 0:
        raise ValueError("License header must be non-empty!")  # noqa: TRY003, EM101
    safe_license_header_lines: list[str] = []
    for line in license_header_contents.split("\n"):
        this_line = f"# {line}" if not line.startswith("#") else line
        safe_license_header_lines.append(this_line)
    return "\n".join(safe_license_header_lines)


def get_staged_files() -> list[str]:
    """Returns list of git staged files."""
    try:
        result = subprocess.run(  # noqa: S603
            ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],  # noqa: S607
            capture_output=True,
            text=True,
            check=True,  # Raise exception on non-zero exit code
        )
    except subprocess.CalledProcessError as e:
        error_msg = f"""Git command failed with code {e.returncode}
        Command: {" ".join(e.cmd)}
        Error output: {e.stderr.strip() or "<no output>"}
        """
        raise RuntimeError(error_msg) from None
    except FileNotFoundError:
        raise RuntimeError("Git executable not found - is Git installed?") from None  # noqa: TRY003, EM101

    if not result.stdout.strip():
        return []
    return result.stdout.splitlines()


@click.command(help="Check that Python files start with a license header.")
@click.option(
    "--check",
    "-c",
    multiple=True,
    type=str,
    help="Either a file or directory. If a directory, then all files that are "
    "accessible will be included (directories)"
    " are searched recursively). Only files that end with *.py will be included. "
    "Acceptable to use multiple "
    "times. All --check files will be included. Must specify at least one *.py file.",
)
@click.option(
    "--modify",
    "-m",
    is_flag=True,
    help="If present, modifies files that don't have the license header. "
    "Otherwise, will error-out if it finds any non-compliant files.",
)
@click.option(
    "--license-header",
    "-l",
    required=False,
    help="If present, loads the license header from this file. Defaults to use " "standard license header.",
)
@click.option(
    "--add-leading",
    "-a",
    is_flag=True,
    help="If present, will ensure that each line of the license header starts with "
    "'#'. "
    "If any line doesn't, then this option will make the program append '# ' to the "
    "start of each line.",
)
@click.option(
    "--replace",
    "-r",
    is_flag=True,
    help="If present, will replace an existing license header. By default, this "
    "program simply appends the"
    "license header text to each .py file. This option will employ a heuristic to "
    "detect if a .py file "
    "starts with a license header. If detected, then this text is removed before the"
    " normal license header appending logic runs.",
)
@click.option(
    "--verbose",
    is_flag=True,
    help="If present, will perform extra (verbose) logging. ",
)
def entrypoint(  # noqa: PLR0915, PLR0912, C901
    check: tuple[str, ...],
    modify: bool,
    license_header: str | None,
    add_leading: bool,
    replace: bool,
    verbose: bool,
) -> None:
    logger.info(f"Files/directories for finding .py files: {check}")  # noqa: G004
    logger.info(f"Modify .py files with license header?:   {modify}")  # noqa: G004
    logger.info(f"Overriding standard license header?:     {license_header}")  # noqa: G004
    logger.info(f"Force each line to start with '# '?:     {add_leading}")  # noqa: G004
    logger.info(f"Check for and replace existing header?:  {replace}")  # noqa: G004
    logger.info(f"Verbose (extra) logging?:                {verbose}")  # noqa: G004
    logger.info("-" * 100)

    if len(check) == 0:
        files_staged = get_staged_files()
        py_staged = []
        for f in files_staged:
            p = Path(f).absolute()
            if p.name.endswith(".py"):
                py_staged.append(f)
        check = tuple(f for f in py_staged)
        logger.info(f"Run license_check.py on the staged files: {check}")  # noqa: G004

    if replace and not modify:
        raise ValueError("Must use --modify if also using --replace !")  # noqa: TRY003, EM101

    if modify and not replace:
        logger.info(
            "WARNING: existing license headers are ignored. "
            "To replace any existing header text, re-run with --replace. "
        )

    # get all files / directories from --check
    files: list[Path] = []
    directories: list[Path] = []
    unknown: list[Path] = []
    for f in check:
        p = Path(f).absolute()
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            directories.append(p)
        else:
            unknown.append(p)

    # check that they all exist
    if len(unknown) > 0:
        raise ValueError(f"Found {len(unknown)} --check things that do not exist!\n\n".join(map(str, unknown)))

    # check that files passed in explicitly from --check end with .py
    non_py_files = [f for f in files if not f.name.endswith(".py")]
    if len(non_py_files) > 0:
        raise ValueError(
            f"Found {len(non_py_files)} files from --check that are not Python files! "
            "(They don't end with .py):\n"
            "\n".join(map(str, non_py_files))
        )

    # resolve the license header: either the default or user override .txt file
    if license_header is not None:
        lic_file = Path(license_header).absolute()
        if not lic_file.is_file():
            raise ValueError(  # noqa: TRY003
                f"Supplied a --license-header, but {license_header} is not a file!"  # noqa: EM102
            )
        with open(str(lic_file)) as rt:  # noqa: PTH123
            license_header_contents: str = rt.read()
        msg_license: str = "Using custom license header"
    else:
        license_header_contents = LICENSE_HEADER
        msg_license = "Using default license header"

    if add_leading:
        msg_license += ". Ensuring each line of the license header starts with '#'"
        license_header_contents = ensure_license_starts_with_pound(license_header_contents)

    logger.info(f"{msg_license}:\n{license_header_contents}" if verbose else f"{msg_license}.")

    # run license check
    try:
        checked_n_files: int = main(
            modify,
            license_header_contents,
            files=files,
            directories=directories,
            replace=replace,
        )
    except ValueError as error:
        logger.info(str(error))
        sys.exit(1)
    else:
        logger.info(
            f"Success! All {checked_n_files} checked have the required license header!"  # noqa: G004
        )


def main(
    modify: bool,
    license_header_contents: str,
    *,
    files: list[Path],
    directories: list[Path],
    replace: bool,
) -> int:
    """Runs license check on all files & files accessible from the directories.

    On failure, raises an error with all noncompliant files. Returns nothing on success.
      See :func:`check_license_project_files` for details.

    Returns the number of files checked on success. On failure, :raises:`ValueError`
      with message containing the # of non-compliant files & their specific
      :class:`LicenseCheckError` errors.

    The :param:`replace` option will heuristically check for an existing license header.
      It will remove and replace this with the :param:`license_header_contents`.
    """
    if len(files) == 0 and len(directories) == 0:
        return 0
    if len(license_header_contents) == 0:
        raise ValueError("Must supply non-empty license header!")  # noqa: TRY003, EM101
    checked = _main(modify, license_header_contents, files, directories, replace)
    if len(checked.noncompliant_files) > 0:
        raise _error(checked.noncompliant_files, checked.n_files, modify)
    return checked.n_files


def _main(
    modify: bool,
    license_header_contents: str,
    files: list[Path],
    directories: list[Path],
    replace: bool,
) -> Checked:
    check_file: Callable[[Path], LicenseCheckError | None] = partial(
        license_check,
        license_header=license_header_contents,
        modify=modify,
        replace=replace,
    )
    check_dir: Callable[[Path], Checked] = partial(
        check_license_project_files,
        modify=modify,
        license_header=license_header_contents,
        replace=replace,
    )

    n_files_checked: int = 0
    noncompliant_files: dict[Path, LicenseCheckError] = {}

    # license check all individual files
    for f in files:
        maybe_err = check_file(f)
        if maybe_err is not None:
            noncompliant_files[f] = maybe_err
    n_files_checked += len(files)

    # license check all directories and their contents, recursively
    for d in directories:
        checked = check_dir(d)
        noncompliant_files.update(checked.noncompliant_files)
        n_files_checked += checked.n_files

    return Checked(noncompliant_files=noncompliant_files, n_files=n_files_checked)


def _error(
    noncompliant_files: Mapping[Path, LicenseCheckError],
    n_files_checked: int,
    modify: bool,
) -> ValueError:
    maybe_modify_msg: str = (
        " You can re-run with '--modify' to automatically add the required" " license header." if not modify else ""
    )
    error_message: str = (
        f"ERROR: There are {len(noncompliant_files)} / {n_files_checked} "
        f"files that do not have the license header!{maybe_modify_msg}\n"
    )
    for pyfile, error in noncompliant_files.items():
        error_message += f"  {pyfile!s}: {error}\n"
    return ValueError(error_message)


if __name__ == "__main__":
    entrypoint()
