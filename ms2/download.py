"""Fetch the pinned public input and reject any change in its bytes."""

import argparse
import hashlib
import os
from pathlib import Path
import tempfile
from urllib.request import urlopen

from .protocol import INPUT_BYTES, INPUT_PATH, SHA256, URL


def verify(path):
    """Check both the length and SHA-256 of an existing input file."""
    path = Path(path)
    if path.stat().st_size != INPUT_BYTES:
        raise ValueError(f"Wrong input length: {path.stat().st_size}; expected {INPUT_BYTES}")
    with path.open("rb") as handle:
        actual = hashlib.file_digest(handle, "sha256").hexdigest()
    if actual != SHA256:
        raise ValueError(f"SHA-256 mismatch: {actual}; expected {SHA256}")
    return actual


def download(path):
    """Write a verified download atomically, or verify the existing file."""
    path = Path(path)
    if path.exists():
        verify(path)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".massspecgym-", dir=path.parent)
    temporary = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as output, urlopen(URL, timeout=60) as response:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
        verify(temporary)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=INPUT_PATH, help="Destination TSV file")
    args = parser.parse_args()
    download(args.output)
    print(f"Verified {args.output}\nSHA-256 {SHA256}")


if __name__ == "__main__":
    main()
