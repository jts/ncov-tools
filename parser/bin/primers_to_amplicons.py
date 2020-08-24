#!/usr/bin/env python
'''
Convert the nCoV primer scheme to a unique amplicon BED file.
'''

import sys
import argparse
import ncov.parser.primers as pr

parser = argparse.ArgumentParser(description='Create amplicon BED file')
parser.add_argument('-p', '--primers', help='Primer scheme in BED format')
parser.add_argument('--offset', default=0, help='Primer offset for coordinates')
parser.add_argument('-o', '--output', default='out.bed',
                    help='filename to write BED to')
parser.add_argument('--bed_type', default='unique_amplicons',
                    help='type of BED to produce (e.g. full, no_primers, unique-amplicons')
if len(sys.argv) <= 1:
    parser.print_help(sys.stderr)
    sys.exit('Invalid number of arguments')
args = parser.parse_args()

primers = pr.read_bed_file(args.primers)
primer_pairs = pr.create_primer_pairs(primers=primers)
amplicon_ranges = pr.create_amplicons(primer_pairs=primer_pairs,
                                      offset=args.offset,
                                      type=args.bed_type)

with open(args.output, 'w') as file_o:
    for line in amplicon_ranges:
        file_o.write('\t'.join(line))
        file_o.write('\n')
file_o.close()
