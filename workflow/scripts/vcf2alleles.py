#!/usr/bin/env python


from ncov.parser.Vcf import *
import argparse
import os
import sys
import re
from copy import deepcopy

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
    assert upper_limit > lower_limit, "VAF upper_limit must be greater than lower_limit"
    if vaf >= upper_limit:
        return False
    elif vaf < upper_limit and vaf >= lower_limit:
        return True
    else:
        sys.exit('Invalid VAF or VAF less than 0.25')
        

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


def separate_mutations(var):
    """
    Several mutations were identified as having multiple
    nucleotide polymorphisms, however the mixture_check.py
    script requires each mutation to be a separate line
    in the alleles.tsv file.  This function returns
    a list of variants.
    """
    var_orig = var
    vars = list()
    tmp_ref = str(var.REF)
    tmp_alt = str(var.ALT[0])
    tmp_pos = int(var.POS)
    for idx, (ref, alt) in enumerate(zip(tmp_ref, tmp_alt)):
        tmp_var = deepcopy(var)
        tmp_var.REF = ref
        tmp_var.ALT[0] = alt
        if ref == alt:
            continue
        else:
            tmp_var.POS = tmp_pos + int(idx)
            vars.append(tmp_var)
    return vars


def is_mnp(var):
    """
    Check if the mutation is a multiple nucleotide polymorphism.
    """
    if len(str(var.ALT[0])) == len(str(var.REF)) and len(str(var.ALT[0])) > 1:
        return True
    elif len(str(var.ALT[0])) != len(str(var.REF)):
        return False


def get_variants_from_vcf(file):
    """
    Returns a list of single nucleotide variants from a VCF file.
    """
    _var_list = list()
    _vcf_reader = vcf.Reader(filename=file)
    for _var in _vcf_reader:
        if is_mnp(var=_var):
            _var_list.extend(separate_mutations(var=_var))
        elif len(str(_var.REF)) == 1 and len(str(_var.ALT[0])) == 1:
            _var_list.append(_var)
        else:
            continue
    return _var_list



def main():
    """
    Main program
    """
    args = init_args()
    var_dict = dict()
    # loop through all VCF files in the path provided
    for file in os.listdir(args.path):
        vcf_list = list()
        if file.endswith(args.pattern):
            vcf_reader = vcf.Reader(filename='/'.join([args.path, file]))
            # import the variants in the VCF file while converting MNPs to
            # SNPs
            for var in vcf_reader:
                if len(str(var.REF)) == 1 and len(str(var.ALT[0])) == 1:
                    vcf_list.append(var)
                # skip indels but process multiple nucleotide polymorphisms
                elif is_mnp(var=var):
                    vcf_list.extend(separate_mutations(var=var))
                else:
                    continue
            # create the alleles dictionary
            for _var in vcf_list: 
                vaf = _var.INFO['VAF'][0]
                if is_ambiguous(vaf):
                    _var.ALT[0] = get_iupac(str(_var.REF), str(_var.ALT[0]))
                var_id = create_vcf_variant_id(var=_var)
                if var_id not in var_dict:
                    var_dict[var_id] = create_default_sample(
                            name=get_sample_name(file=file, pattern=args.pattern),
                            pos=_var.POS,
                            ref=_var.REF,
                            alt=_var.ALT[0])
                elif var_id in var_dict:
                    var_dict[var_id]['samples'].append(
                            {'name': get_sample_name(file=file, pattern=args.pattern),
                             'pos': _var.POS, 'ref': _var.REF, 'alt': _var.ALT[0]})
    print('\t'.join(['name', 'pos', 'ref_allele', 'alt_allele', 'samples_with_allele']))
    for item in var_dict:
        write_alleles_line(var_dict[item])


if __name__ == '__main__':
    main()


#__END__

