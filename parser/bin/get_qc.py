#!/usr/bin/env python
'''
A Python package for summarizing QC data from the ncov-tools pipeline.
'''

import argparse
import sys
import ncov.parser.qc as qc
import ncov.parser

parser = argparse.ArgumentParser(description="Tool for summarizing QC data")
parser.add_argument('-c', '--consensus', help='<sample>.consensus.fasta file to process')
parser.add_argument('-v', '--variants',
                    help='<sample>.variants.tsv file to process')
parser.add_argument('-e', '--coverage',
                    help='<sample>.per_base_coverage.bed file to process')
parser.add_argument('-i', '--indel', action='store_true',
                    help='flag to determine whether to count indels')
parser.add_argument('-m', '--meta', default=None,
                    help='full path to the metadata file')
parser.add_argument('-a', '--alleles',
                    help='full path to the alleles.tsv file')
parser.add_argument('-s', '--sample',
                    help='name of sample being processed')
parser.add_argument('-p', '--platform', default='illumina',
                    help='sequencing platform used')
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit('Invalid number of arguments')
args = parser.parse_args()

qc_line = dict()
qc_line.update({'sample' : args.sample})

try:
    meta = ncov.parser.Meta(file=args.meta)
    meta.import_metadata()
    qc_line.update(meta.data[args.sample])
except:
    qc_line.update({'qpcr_ct' : 'NA', 'collection_date' : 'NA',
                    'num_months' : 'NA', 'num_weeks' : 'NA'})

if args.platform == 'illumina':
    if str(args.variants).endswith('.variants.tsv'):
        vars = ncov.parser.Variants(file=args.variants)
        qc_line.update(vars.get_total_variants())
    else:
        sys.exit('Must be a valid variant.tsv file for the Illumina platform')
elif args.platform == 'oxford-nanopore':
    if str(args.variants).endswith('.vcf'):
        vars = ncov.parser.Vcf(file=args.variants)
        qc_line.update(vars.get_variant_counts())
    else:
        sys.exit('Must be a valid VCF file for the Oxford-Nanopore platform')

alleles = ncov.parser.Alleles(file=args.alleles)
qc_line.update(alleles.get_variant_counts(sample=args.sample))

cons = ncov.parser.Consensus(file=args.consensus)
qc_line.update(cons.count_iupac_in_fasta())
qc_line.update(cons.get_genome_completeness())

coverage = ncov.parser.PerBaseCoverage(file=args.coverage)
qc_line.update(coverage.get_coverage_stats())
qc_line.update({'qc_pass' : 'NA'})

qc.write_qc_summary_header()
qc.write_qc_summary(summary=qc_line)
