#!/usr/bin/env python3
"""Filter FASTA sequences to standard amino acids (no U/O/X/etc).

Produces two files next to the input:
 - <base>_standard_only.fa : sequences composed entirely of the 20 standard AAs
 - <base>_standard_cleaned.fa : sequences with non-standard residues removed (skips empty results)
"""
import os
import argparse

STANDARD = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy')

def filter_standard(input_fasta, out_only=None, out_cleaned=None):
    input_fasta = os.path.abspath(os.path.expanduser(input_fasta))
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    base = os.path.splitext(os.path.basename(input_fasta))[0]
    dirn = os.path.dirname(input_fasta)
    if out_only is None:
        out_only = os.path.join(dirn, f"{base}_standard_only.fa")
    if out_cleaned is None:
        out_cleaned = os.path.join(dirn, f"{base}_standard_cleaned.fa")

    total = kept_only = cleaned_written = removed_nonstandard = 0

    with open(input_fasta, 'r') as fin, open(out_only, 'w') as fout_only, open(out_cleaned, 'w') as fout_clean:
        header = None
        seq_lines = []
        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    total += 1
                    seq = ''.join(seq_lines)
                    if seq and all((c in STANDARD) for c in seq):
                        fout_only.write(header + '\n')
                        for i in range(0, len(seq), 60):
                            fout_only.write(seq[i:i+60] + '\n')
                        kept_only += 1
                    # cleaned: remove non-standard
                    cleaned = ''.join([c for c in seq if c in STANDARD])
                    if cleaned:
                        fout_clean.write(header + '\n')
                        for i in range(0, len(cleaned), 60):
                            fout_clean.write(cleaned[i:i+60] + '\n')
                        cleaned_written += 1
                    else:
                        removed_nonstandard += 1

                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # last record
        if header is not None:
            total += 1
            seq = ''.join(seq_lines)
            if seq and all((c in STANDARD) for c in seq):
                fout_only.write(header + '\n')
                for i in range(0, len(seq), 60):
                    fout_only.write(seq[i:i+60] + '\n')
                kept_only += 1
            cleaned = ''.join([c for c in seq if c in STANDARD])
            if cleaned:
                fout_clean.write(header + '\n')
                for i in range(0, len(cleaned), 60):
                    fout_clean.write(cleaned[i:i+60] + '\n')
                cleaned_written += 1
            else:
                removed_nonstandard += 1

    return {
        'input': input_fasta,
        'out_only': out_only,
        'out_cleaned': out_cleaned,
        'total': total,
        'kept_only': kept_only,
        'cleaned_written': cleaned_written,
        'removed_nonstandard': removed_nonstandard,
    }

def main():
    p = argparse.ArgumentParser(description='Filter FASTA for standard amino acids')
    p.add_argument('input', nargs='?', default=os.path.join(os.path.dirname(__file__), '..', 'hs_cytoplasmic', 'hs_cytoplasm_filtered_lt5000.fa'))
    args = p.parse_args()

    res = filter_standard(args.input)
    print(f"Processed: {res['input']}")
    print(f"Total sequences: {res['total']}")
    print(f"Sequences fully standard (written to): {res['kept_only']} -> {res['out_only']}")
    print(f"Sequences cleaned (non-standard removed) written: {res['cleaned_written']} -> {res['out_cleaned']}")
    print(f"Sequences removed entirely due to non-standard residues: {res['removed_nonstandard']}")

if __name__ == '__main__':
    main()
