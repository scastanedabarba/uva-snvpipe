#!/usr/bin/env python3
"""Basic FASTQ read statistics; output format differs from legacy Rivanna fastqstats."""
import gzip
import sys

def main(path):
    opener = gzip.open if path.endswith('.gz') else open
    reads = bases = 0
    with opener(path, 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not seq or not plus or not qual or not header.startswith('@') or not plus.startswith('+') or len(seq.strip()) != len(qual.strip()):
                raise ValueError(f'Malformed FASTQ record at read {reads+1}')
            reads += 1
            bases += len(seq.strip())
    print('metric\tvalue')
    print(f'reads\t{reads}')
    print(f'bases\t{bases}')
    print(f'mean_read_length\t{bases/reads:.2f}' if reads else 'mean_read_length\t0')

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('Usage: fastqstats.py input.fastq[.gz]')
    main(sys.argv[1])
