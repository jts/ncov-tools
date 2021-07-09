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


def is_ambiguous(vaf, upper_limit=0.75, lower_limit=0.25):
    """
    Determine whether the ALT allele should be an IUPAC code
    based on the VAF.
    """
    assert vaf > lower_limit, "VAF must be greater than 0.25"
    if vaf >= upper_limit:
        return False
    elif vaf < upper_limit:
        return True
    else:
        sys.exit('Invalid VAF')
        

def get_iupac(ref, alt):
    """
    From the concatenated reference and alternate allele,
    return the corresponding IUPAC code.
    """
    ref_alt = ''.join([ref, alt])
    iupac_map = { "AC" : "M",
                  "CA" : "M",
                  "AG" : "R",
                  "GA" : "R",
                  "AT" : "W",
                  "TA" : "W",
                  "CG" : "S",
                  "GC" : "S",
                  "CT" : "Y",
                  "TC" : "Y",
                  "TG" : "K",
                  "GT" : "K" }
    return iupac_map[ref_alt]


def main():
    args = init_args()
    var_dict = dict()
    for file in os.listdir(args.path):
        if file.endswith(args.pattern):
            vcf_reader = vcf.Reader(filename='/'.join([args.path, file]))
            for var in vcf_reader:
                # skip indels
                if len(str(var.REF)) > 1 or len(str(var.ALT[0])) > 1:
                    continue
                vaf = var.INFO['VAF'][0]
                if is_ambiguous(vaf):
                    var.ALT[0] = get_iupac(str(var.REF), str(var.ALT[0]))
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

