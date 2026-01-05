#!/usr/bin/env python3
"""Extract FASTA entries whose accession lacks a matching results subdirectory.

Usage:
  python scripts/find_missing_accessions.py input.fa /path/to/output_dirs [-o missing.fa]

The directory is expected to contain subdirectories named like:
  00012345_Q8WZ33
The accession is taken from the portion after the first underscore.
FASTA headers are assumed to start with `>` followed by the accession token.
"""
import argparse
import os
import sys
from typing import Iterable, List, Set, Tuple

Entry = Tuple[str, str]


def read_fasta(path: str) -> List[Entry]:
    """Return list of (header, sequence) entries from a FASTA file."""
    entries: List[Entry] = []
    header: str | None = None
    seq_lines: List[str] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            entries.append((header, "".join(seq_lines)))
    return entries


def write_fasta(path: str | None, entries: Iterable[Entry]) -> None:
    """Write entries to a file or stdout."""
    out = sys.stdout if path is None else open(path, "w")
    try:
        for header, seq in entries:
            out.write(header + "\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")
    finally:
        if path is not None:
            out.close()


def accession_from_header(header: str) -> str:
    """Extract accession token from a FASTA header (after '>')."""
    return header[1:].split()[0]


def accessions_from_dirs(base_dir: str) -> Set[str]:
    """Collect accessions from subdirectories named '<number>_<accession>'."""
    accs: Set[str] = set()
    for entry in os.scandir(base_dir):
        if not entry.is_dir():
            continue
        name = entry.name
        if "_" not in name:
            continue
        _, acc = name.split("_", 1)
        if acc:
            accs.add(acc)
    return accs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Emit FASTA entries whose accession lacks a corresponding subdirectory."
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument(
        "directory", help="Directory containing subdirectories named like '12345_ACCESSION'"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA path (defaults to stdout)",
        default=None,
    )
    args = parser.parse_args()

    entries = read_fasta(args.fasta)
    dir_accessions = accessions_from_dirs(args.directory)

    missing = []
    for header, seq in entries:
        acc = accession_from_header(header)
        if acc not in dir_accessions:
            missing.append((header, seq))

    write_fasta(args.output, missing)

    print(
        f"Total sequences: {len(entries)} | Missing directories: {len(missing)}",
        file=sys.stderr,
    )
    if args.output:
        print(f"Wrote missing sequences to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
