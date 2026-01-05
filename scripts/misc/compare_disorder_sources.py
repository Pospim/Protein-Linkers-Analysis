#!/usr/bin/env python3
"""
Compare AUIPROD disorder predictions against NetSurfP disorder outputs.

By default this script looks for:
- AUIPROD scores in hs_cytoplasmic/disorder.json
- NetSurfP outputs in hs_cytoplasmic/netsurfp/<run>_<accession>/<run>_<accession>.json

It reports per-accession correlation and error metrics and can optionally write
the per-accession table to CSV.
"""

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Optional, Tuple


def load_auiprod(path: Path) -> Dict[str, List[float]]:
    """Return {accession: [disorder scores]} from the AUIPROD JSON file."""
    with path.open() as f:
        data = json.load(f)
    return {acc: [float(v) for v in vals] for acc, vals in data.items()}


def load_netsurfp_disorder(
    root: Path, limit: Optional[int] = None
) -> Dict[str, List[float]]:
    """
    Return {accession: [disorder scores]} from NetSurfP output folders.

    Each NetSurfP folder is expected to contain a single .json file with a
    "disorder" array. The accession is taken from the JSON "desc" field when
    present, otherwise from the folder name suffix after the first underscore.
    """
    disorders: Dict[str, List[float]] = {}
    entries = sorted(p for p in root.iterdir() if p.is_dir())
    if limit is not None:
        entries = entries[:limit]

    for entry in entries:
        json_files = sorted(entry.glob("*.json"))
        if not json_files:
            continue

        with json_files[0].open() as f:
            payload = json.load(f)

        accession = (payload.get("desc") or entry.name.split("_", 1)[-1]).strip()
        disorder_vals = payload.get("disorder")
        if disorder_vals is None:
            continue
        disorders[accession] = [float(v) for v in disorder_vals]
    return disorders


def pearson(x: List[float], y: List[float]) -> float:
    """Compute the Pearson correlation between two equal-length lists."""
    n = len(x)
    if n == 0 or n != len(y):
        return math.nan
    mx = mean(x)
    my = mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = sum((a - mx) ** 2 for a in x)
    dy = sum((b - my) ** 2 for b in y)
    if dx <= 0 or dy <= 0:
        return math.nan
    return num / math.sqrt(dx * dy)


def compare(
    auiprod: Dict[str, List[float]],
    netsurfp: Dict[str, List[float]],
) -> Tuple[
    List[Dict[str, float]],
    List[str],
    List[str],
    List[Tuple[str, int, int]],
]:
    """
    Compare AUIPROD and NetSurfP disorder scores.

    Returns:
    - per-accession metrics (list of dicts)
    - accessions missing NetSurfP outputs
    - accessions missing AUIPROD outputs
    - length mismatches (accession, auiprod_len, netsurfp_len)
    """
    results: List[Dict[str, float]] = []
    missing_netsurfp: List[str] = []
    length_mismatches: List[Tuple[str, int, int]] = []

    for acc, au_scores in auiprod.items():
        net_scores = netsurfp.get(acc)
        if net_scores is None:
            missing_netsurfp.append(acc)
            continue

        if len(au_scores) != len(net_scores):
            length_mismatches.append((acc, len(au_scores), len(net_scores)))
            continue

        p = pearson(au_scores, net_scores)
        mae = sum(abs(a - b) for a, b in zip(au_scores, net_scores)) / len(au_scores)
        rmse = math.sqrt(
            sum((a - b) ** 2 for a, b in zip(au_scores, net_scores)) / len(au_scores)
        )

        results.append(
            {
                "accession": acc,
                "length": len(au_scores),
                "pearson": p,
                "mae": mae,
                "rmse": rmse,
                "auiprod_mean": mean(au_scores),
                "netsurfp_mean": mean(net_scores),
            }
        )

    missing_auiprod = sorted(set(netsurfp.keys()) - set(auiprod.keys()))
    return results, missing_netsurfp, missing_auiprod, length_mismatches


def write_csv(rows: Iterable[Dict[str, float]], path: Path) -> None:
    fieldnames = [
        "accession",
        "length",
        "pearson",
        "mae",
        "rmse",
        "auiprod_mean",
        "netsurfp_mean",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare AUIPROD disorder predictions to NetSurfP outputs."
    )
    parser.add_argument(
        "--auiprod",
        default=Path("hs_cytoplasmic") / "disorder.json",
        type=Path,
        help="Path to AUIPROD disorder JSON file.",
    )
    parser.add_argument(
        "--netsurfp",
        default=Path("hs_cytoplasmic") / "netsurfp",
        type=Path,
        help="Path to NetSurfP output directory.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional path to write per-accession metrics as CSV.",
    )
    parser.add_argument(
        "--show-top",
        type=int,
        default=5,
        help="Show the top N accessions by Pearson correlation (default: 5).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on how many NetSurfP folders to process (for quick tests).",
    )
    args = parser.parse_args()

    auiprod = load_auiprod(args.auiprod)
    netsurfp = load_netsurfp_disorder(args.netsurfp, limit=args.limit)
    results, missing_netsurfp, missing_auiprod, length_mismatches = compare(
        auiprod, netsurfp
    )

    print(
        f"NetSurfP accessions loaded: {len(netsurfp)} "
        f"(limit={'all' if args.limit is None else args.limit})"
    )
    pearsons = [r["pearson"] for r in results if not math.isnan(r["pearson"])]
    maes = [r["mae"] for r in results]

    missing_note = (
        " (computed against limited NetSurfP subset)" if args.limit else ""
    )

    print(f"Matched accessions: {len(results)}")
    print(f"Missing NetSurfP files: {len(missing_netsurfp)}{missing_note}")
    print(f"Missing AUIPROD records: {len(missing_auiprod)}")
    print(f"Length mismatches: {len(length_mismatches)}")

    if pearsons:
        print(
            f"Pearson r: mean={mean(pearsons):.4f}, "
            f"median={median(pearsons):.4f}, "
            f"min={min(pearsons):.4f}, max={max(pearsons):.4f}, n={len(pearsons)}"
        )
    if maes:
        print(f"MAE: mean={mean(maes):.4f}, median={median(maes):.4f}")

    if args.show_top and results:
        top = sorted(
            [r for r in results if not math.isnan(r["pearson"])],
            key=lambda r: r["pearson"],
            reverse=True,
        )[: args.show_top]
        print("\nTop correlations:")
        for r in top:
            print(
                f"  {r['accession']}: r={r['pearson']:.4f}, "
                f"mae={r['mae']:.4f}, rmse={r['rmse']:.4f}, len={r['length']}"
            )

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        write_csv(results, args.csv)
        print(f"\nPer-accession metrics written to: {args.csv}")


if __name__ == "__main__":
    main()
