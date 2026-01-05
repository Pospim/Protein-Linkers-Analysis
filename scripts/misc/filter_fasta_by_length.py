#!/usr/bin/env python3
"""Filter a FASTA file by sequence length.

Writes sequences with length strictly less than `max_length` to output FASTA.
"""
import argparse
import os

def filter_fasta(input_fasta, output_fasta, max_length=5000):
    input_fasta = os.path.abspath(os.path.expanduser(input_fasta))
    output_fasta = os.path.abspath(os.path.expanduser(output_fasta))

    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    kept = 0
    total = 0

    with open(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
        header = None
        seq_lines = []
        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    total += 1
                    seq = ''.join(seq_lines)
                    if len(seq) < max_length:
                        fout.write(header + '\n')
                        # wrap at 60 chars
                        for i in range(0, len(seq), 60):
                            fout.write(seq[i:i+60] + '\n')
                        kept += 1
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # last record
        if header is not None:
            total += 1
            seq = ''.join(seq_lines)
            if len(seq) < max_length:
                fout.write(header + '\n')
                for i in range(0, len(seq), 60):
                    fout.write(seq[i:i+60] + '\n')
                kept += 1

    print(f"Filtered FASTA written to: {output_fasta}")
    print(f"Total sequences: {total}, kept (len<{max_length}): {kept}, removed: {total - kept}")


def main():
    parser = argparse.ArgumentParser(description='Filter FASTA by max sequence length')
    parser.add_argument('input', nargs='?', default=os.path.join(os.path.dirname(__file__), '..', 'hs_cytoplasmic', 'hs_cytoplasmic.fa'))
    parser.add_argument('-o', '--output', help='Output FASTA path')
    parser.add_argument('-m', '--max-length', type=int, default=5000, help='Maximum sequence length (exclusive)')
    args = parser.parse_args()

    input_fasta = args.input
    out = args.output
    if not out:
        inp_dir = os.path.dirname(os.path.abspath(os.path.expanduser(input_fasta)))
        base = os.path.splitext(os.path.basename(input_fasta))[0]
        out = os.path.join(inp_dir, f"{base}_filtered_lt{args.max_length}.fa")

    filter_fasta(input_fasta, out, max_length=args.max_length)

if __name__ == '__main__':
    main()
