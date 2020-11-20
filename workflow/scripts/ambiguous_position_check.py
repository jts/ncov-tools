#!/usr/bin/env python
'''
Parse an alleles.tsv file to identify possible problematic positions
'''

import argparse
import sys
import ncov.parser
from collections import defaultdict

parser = argparse.ArgumentParser(description="Generate a report on genomic positions that have frequent ambiguous variants")
parser.add_argument('-a', '--alleles',
                    help='full path to the alleles.tsv file')
parser.add_argument('-m', '--min-count', type=int, default=3,
                    help='Report all positions that have more than m ambiguous observations')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit('Invalid number of arguments')
args = parser.parse_args()

alleles = ncov.parser.Alleles(file=args.alleles)

ambiguous_counts = defaultdict(int)
ambiguous_alleles = defaultdict(dict)
for sample, records in alleles.data.items():
    for position in records:
        aa = records[position]['alt']
        if aa != 'N' and ncov.parser.is_variant_iupac(aa):
            p = int(position)
            ambiguous_counts[p] += 1
            ambiguous_alleles[p][aa] = 1

print("position\tcount\talleles")
for position in sorted(ambiguous_counts.keys()):
    count = ambiguous_counts[position]
    if count >= args.min_count:
        print("%d\t%d\t%s" % (int(position), count, ",".join(ambiguous_alleles[position])))
