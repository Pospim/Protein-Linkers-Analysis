#!/usr/bin/env python3
"""Split a FASTA into 4 batches with ~equal numbers of sequences.
Writes files to the same `hs_cytoplasmic/` folder as `batch_1.fa`..`batch_4.fa`.
"""
import os

ROOT = os.path.dirname(os.path.dirname(__file__))
INPUT = os.path.join(ROOT, 'hs_cytoplasmic', 'hs_cytoplasm_filtered_lt5000_standard_only.fa')
OUT_DIR = os.path.join(ROOT, 'hs_cytoplasmic')
OUT_FILES = [os.path.join(OUT_DIR, f'batch_{i+1}.fa') for i in range(4)]

def read_fasta(path):
    entries = []
    header = None
    seq_lines = []
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    entries.append((header, ''.join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            entries.append((header, ''.join(seq_lines)))
    return entries


def write_batch(path, entries):
    with open(path, 'w') as fh:
        for header, seq in entries:
            fh.write(header + '\n')
            # write sequence in lines of up to 80 chars for readability
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')


def main():
    entries = read_fasta(INPUT)
    n = len(entries)
    if n == 0:
        print('No sequences found in', INPUT)
        return
    base = n // 4
    rem = n % 4
    sizes = [base + (1 if i < rem else 0) for i in range(4)]

    idx = 0
    created = []
    for batch_idx, size in enumerate(sizes):
        batch_entries = entries[idx: idx+size]
        outpath = OUT_FILES[batch_idx]
        write_batch(outpath, batch_entries)
        created.append((outpath, len(batch_entries)))
        idx += size

    print('Input:', INPUT)
    print('Total sequences:', n)
    for p, c in created:
        print(os.path.relpath(p, ROOT) + ':', c)

if __name__ == '__main__':
    main()
