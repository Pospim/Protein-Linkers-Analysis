import csv
import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import make_dssp_dict


# Map DSSP 8-state SS codes to 3-state
def dssp_to_3state(code: str) -> str:
    if code in {"H", "G", "I"}:
        return "H"
    if code in {"E", "B"}:
        return "E"
    return "C"  # coil, turn, bend


def per_residue_plddt(model) -> dict:
    plddt = {}
    for chain in model:
        chain_id = chain.id
        for residue in chain:
            hetflag, resseq, icode = residue.get_id()
            if hetflag != " ":
                continue
            if "CA" not in residue:
                continue
            key = (chain_id, residue.get_id())
            plddt[key] = float(residue["CA"].bfactor)
    return plddt


def run_mkdssp_to_file(pdb_path: str, dssp_out_path: str, mkdssp_exe: str = "mkdssp") -> None:
    """
    Run mkdssp explicitly and force classic DSSP output.
    This avoids Biopython trying to parse/guess mmCIF output formats.
    """
    exe = shutil.which(mkdssp_exe) or mkdssp_exe
    # Force classic DSSP output explicitly
    cmd = [exe, "--output-format", "dssp", pdb_path, dssp_out_path]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)


def process_pdb_gz(pdb_gz_path: Path, out_dir: Path, mkdssp_exe: str = "mkdssp") -> Path:
    with gzip.open(pdb_gz_path, "rt") as fh:
        pdb_text = fh.read()

    tmp_pdb = None
    tmp_dssp = None
    try:
        # Write temporary PDB
        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp:
            tmp.write(pdb_text)
            tmp_pdb = tmp.name

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_gz_path.stem, tmp_pdb)
        model = structure[0]

        # 1) pLDDT per residue (Cα)
        plddt_map = per_residue_plddt(model)

        # 2) Run mkdssp ourselves -> classic DSSP file, then parse via Biopython helper
        with tempfile.NamedTemporaryFile("w", suffix=".dssp", delete=False) as tmp:
            tmp_dssp = tmp.name

        run_mkdssp_to_file(tmp_pdb, tmp_dssp, mkdssp_exe=mkdssp_exe)

        # Parse DSSP output file to a dictionary
        # Keys look like: (chain_id, (' ', resseq, icode))
        dssp_dict = make_dssp_dict(tmp_dssp)[0]

        out_path = out_dir / f"{pdb_gz_path.name}.ss_plddt.csv"
        with open(out_path, "w", newline="") as out_f:
            w = csv.writer(out_f)
            w.writerow([
                "file",
                "chain",
                "resseq",
                "icode",
                "aa",
                "dssp8",
                "ss3",
                "plddt_ca",
                "disorder_plddt",
            ])

            for (chain_id, res_id), vals in dssp_dict.items():
                _, resseq, icode = res_id

                # make_dssp_dict returns tuples where:
                # vals[0] = amino acid (one-letter)
                # vals[1] = DSSP secondary structure (8-state)
                aa = vals[0]
                dssp8 = vals[1]
                ss3 = dssp_to_3state(dssp8)

                plddt_val = plddt_map.get((chain_id, res_id))
                dis = ""
                if plddt_val is not None:
                    dis = plddt_val / 100.0

                w.writerow([
                    pdb_gz_path.name,
                    chain_id,
                    int(resseq),
                    (icode.strip() if isinstance(icode, str) else ""),
                    aa,
                    dssp8,
                    ss3,
                    (f"{plddt_val:.2f}" if plddt_val is not None else ""),
                    dis,
                ])

        return out_path

    finally:
        for p in (tmp_pdb, tmp_dssp):
            if p and os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Batch DSSP + per-residue pLDDT (CA) from AlphaFold PDB.GZ files")
    ap.add_argument("--input_dir", required=True, help="Directory containing *.pdb.gz AlphaFold files")
    ap.add_argument("--out", default="ss_plddt_out", help="Output directory")
    ap.add_argument("--mkdssp", default="mkdssp", help="Path to mkdssp executable (default: mkdssp on PATH)")
    args = ap.parse_args()

    inp = Path(args.input_dir)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_gzs = sorted(inp.glob("*.pdb.gz"))
    if not pdb_gzs:
        raise SystemExit(f"No *.pdb.gz files found in: {inp}")

    for p in pdb_gzs:
        out_path = process_pdb_gz(p, out_dir, mkdssp_exe=args.mkdssp)
        print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
