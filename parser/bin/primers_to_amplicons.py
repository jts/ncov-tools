#!/usr/bin/env python
'''
Convert the nCoV primer scheme to a unique amplicon BED file.
'''

import sys
import argparse
import ncov.parser.primers as pr

parser = argparse.ArgumentParser(description='Create unique amplicon BED file')
parser.add_argument('-p', '--primers', help='Primer scheme in BED format')
parser.add_argument('--offset', default=30, help='Primer offset for coordinates')
parser.add_argument('-o', '--output', default='out.bed',
                    help='filename to write BED to')
parser.add_argument('--pattern', default='nCoV-2019_',
                    help='amplicon name pattern')
parser.add_argument('--remove_primers', action='store_true',
                    help='remove primer sequences from amplicons')
if len(sys.argv) <= 1:
    parser.print_help(sys.stderr)
    sys.exit('Invalid number of arguments')
args = parser.parse_args()

primers = pr.read_bed_file(args.primers)
primer_pairs = pr.create_primer_pairs(primers=primers)
amplicon_ranges = pr.create_amplicon_range(primer_pairs=primer_pairs,
                                           pattern=args.pattern)
if args.remove_primers:
    amplicons = pr.create_unique_amplicons(amplicons=amplicon_ranges,
                                                  offset=args.offset)
else:
    amplicons = amplicon_ranges

with open(args.output, 'w') as file_o:
    for line in amplicons:
        file_o.write('\t'.join(line))
        file_o.write('\n')
file_o.close()
