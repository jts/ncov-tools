#!/usr/bin/env python


from ncov.parser.Vcf import *
import argparse
import os
import sys
import re


def init_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', required=True,
                        help='path containing VCF files to process')
    parser.add_argument('-t', '--pattern', default='.variants.norm.vcf',
                        help='the VCF file pattern [default: .variants.norm.vcf]')
    return parser.parse_args()

def create_default_sample(name, pos, ref, alt):
    return {'samples' : [{'name' : name, 'pos': pos, 'ref': ref, 'alt': alt}]}


def get_sample_name(file, pattern):
    return re.sub(pattern, '', file)


def write_alleles_line(row):
    """
    Write out sample mutation as follows:
    * name
    * pos
    * ref_allele
    * alt_allele
    * samples_with_allele
    """
    samples_with_alleles = len(row['samples'])
    for sample in row['samples']:
        if samples_with_alleles < 1:
            continue
        # ignore indels
        elif len(str(sample['ref'])) > 1 or len(str(sample['alt'])) > 1:
            continue
        else:
            print('\t'.join([str(sample['name']),
                             str(sample['pos']),
                             str(sample['ref']),
                             str(sample['alt']),
                             str(samples_with_alleles)]))


def main():
    args = init_args()
    var_dict = dict()
    for file in os.listdir(args.path):
        if file.endswith(args.pattern):
            vcf_reader = vcf.Reader(filename='/'.join([args.path, file]))
            for var in vcf_reader:
                var_id = create_vcf_variant_id(var=var)
                if var_id not in var_dict:
                    var_dict[var_id] = create_default_sample(name=get_sample_name(file=file, pattern=args.pattern), pos=var.POS, ref=var.REF, alt=var.ALT[0])
                elif var_id in var_dict:
                    var_dict[var_id]['samples'].append({'name': get_sample_name(file=file, pattern=args.pattern), "pos": var.POS, "ref": var.REF, "alt": var.ALT[0]})
    print('\t'.join(['name', 'pos', 'ref_allele', 'alt_allele', 'samples_with_allele']))
    for item in var_dict:
        write_alleles_line(var_dict[item])


if __name__ == '__main__':
    main()

