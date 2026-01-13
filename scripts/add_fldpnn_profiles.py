#!/usr/bin/env python3
"""Add flDPnn disorder + PSIPRED SS to per-residue profiles.

This script:
- streams hs_cytoplasmic/fldpnn.json (large JSON array)
- keeps only accessions present in filtered_domains_disorders.json
- computes per-residue disorder (flDPnn_score) and SS3 from PSIPRED helix/strand
- writes updated per_residue_profiles.json and a mean-disorder summary
"""

import argparse
import json
import os
from typing import Dict, Iterable, List, Tuple


def iter_json_array(path: str) -> Iterable[dict]:
    """Stream a top-level JSON array from disk without loading it all."""
    decoder = json.JSONDecoder()
    buf = ""
    in_array = False

    with open(path, "r", encoding="utf-8") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            buf += chunk
            while True:
                if not in_array:
                    idx = buf.find("[")
                    if idx == -1:
                        buf = buf[-1:]  # keep a tiny tail in case '[' splits
                        break
                    buf = buf[idx + 1 :]
                    in_array = True

                buf = buf.lstrip()
                if not buf:
                    break
                if buf[0] == "]":
                    return

                try:
                    obj, idx = decoder.raw_decode(buf)
                except json.JSONDecodeError:
                    break
                yield obj
                buf = buf[idx:].lstrip()
                if buf.startswith(","):
                    buf = buf[1:]


def parse_csv_floats(value, field_name: str, acc: str, debug: bool) -> List[float]:
    if value is None:
        return []
    if isinstance(value, list):
        values = []
        for idx, v in enumerate(value):
            if v is None or str(v).upper() == "NULL":
                if debug:
                    print(f"[warn] {acc} {field_name}[{idx}] is NULL; skipping")
                continue
            try:
                values.append(float(v))
            except ValueError:
                if debug:
                    print(f"[warn] {acc} {field_name}[{idx}]='{v}' not float; skipping")
        return values
    text = str(value).strip()
    if not text:
        return []
    parts = [p for p in text.split(",") if p]
    values = []
    for idx, p in enumerate(parts):
        if p.upper() == "NULL":
            if debug:
                print(f"[warn] {acc} {field_name}[{idx}] is NULL; skipping")
            continue
        try:
            values.append(float(p))
        except ValueError:
            if debug:
                print(f"[warn] {acc} {field_name}[{idx}]='{p}' not float; skipping")
    return values


def normalize_prob(val: float) -> float:
    if val is None:
        return 0.0
    if val > 1.0:
        return val / 100.0
    return val


def derive_ss3(helix_vals: List[float], strand_vals: List[float]) -> List[str]:
    n = min(len(helix_vals), len(strand_vals))
    ss3 = []
    for i in range(n):
        ph = normalize_prob(helix_vals[i])
        pe = normalize_prob(strand_vals[i])
        pc = 1.0 - ph - pe
        if pc < 0:
            pc = 0.0
        if ph >= pe and ph >= pc:
            ss3.append("H")
        elif pe >= ph and pe >= pc:
            ss3.append("E")
        else:
            ss3.append("C")
    return ss3


def mean(values: List[float]):
    if not values:
        return None
    return sum(values) / len(values)


def load_json(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def main() -> None:
    parser = argparse.ArgumentParser(description="Add flDPnn profiles to per-residue data")
    parser.add_argument(
        "--fldpnn",
        default="~/Desktop/work/protein_linkers/hs_cytoplasmic/fldpnn.json",
        help="Path to fldpnn.json (large JSON array)",
    )
    parser.add_argument(
        "--domains",
        default="~/Desktop/work/protein_linkers/hs_cytoplasmic/filtered_domains_disorders.json",
        help="Path to filtered_domains_disorders.json",
    )
    parser.add_argument(
        "--profiles",
        default="~/Desktop/work/protein_linkers/hs_cytoplasmic/per_residue_profiles.json",
        help="Path to per_residue_profiles.json",
    )
    parser.add_argument(
        "--out-profiles",
        default="~/Desktop/work/protein_linkers/hs_cytoplasmic/per_residue_profiles_with_fldpnn.json",
        help="Output path for updated per_residue_profiles",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug messages while parsing",
    )
    args = parser.parse_args()

    fldpnn_path = os.path.expanduser(args.fldpnn)
    domains_path = os.path.expanduser(args.domains)
    profiles_path = os.path.expanduser(args.profiles)
    out_profiles = os.path.expanduser(args.out_profiles)

    domains = load_json(domains_path)
    profiles = load_json(profiles_path)
    domain_accessions = set(domains.keys())

    updated = 0

    for entry in iter_json_array(fldpnn_path):
        acc = entry.get("ACC")
        if not acc or acc not in domain_accessions:
            continue

        if args.debug:
            print(f"[info] processing {acc}")

        dis_vals = parse_csv_floats(entry.get("flDPnn_score"), "flDPnn_score", acc, args.debug)
        helix_vals = parse_csv_floats(entry.get("PSIPRED_helix"), "PSIPRED_helix", acc, args.debug)
        strand_vals = parse_csv_floats(entry.get("PSIPRED_strand"), "PSIPRED_strand", acc, args.debug)
        ss3 = derive_ss3(helix_vals, strand_vals)

        if args.debug:
            print(
                f"[info] {acc} lengths: dis={len(dis_vals)} helix={len(helix_vals)} "
                f"strand={len(strand_vals)} ss3={len(ss3)}"
            )

        if acc not in profiles:
            profiles[acc] = {}

        profiles[acc]["fldpnn_dis"] = dis_vals
        profiles[acc]["fldpnn_ss"] = ss3
        updated += 1

    with open(out_profiles, "w", encoding="utf-8") as f:
        json.dump(profiles, f, indent=2)

    print(f"Updated {updated} accessions from fldpnn")
    print(f"Wrote: {out_profiles}")


if __name__ == "__main__":
    main()
