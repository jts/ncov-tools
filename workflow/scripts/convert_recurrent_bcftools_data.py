#!/usr/bin/env python


import argparse
import os
import sys
import vcf
import re


def init_args():
    """
    Initialize command line arguments
    """
    description = 'Process bcftools vcf and generate a table of recurrent mutations'
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True,
            help='aa table file to convert')
    parser.add_argument('-s', '--sample', required=True,
            help='name of sample to process')
    parser.add_argument('-o', '--output', required=True,
            help='filename for output')
    return parser.parse_args()


def convert_aa(bcsq):
    """
    Convert the BCFTools amino acid change
    to AposB format (similar to SNPEff)
    """
    pass


def get_aa_from_info(var, index=5):
    """
    Get the amino acid (aa) change in the CSQ annotation from
    the VCF INFO string
    """
    aa = str()
    if 'BCSQ' in var.INFO:
        aa = var.INFO['BCSQ'][0]
        return aa.split('|')[index]
    else:
        return None


def get_aa_components(aa):
    """
    The AA string can be broken down into the position
    and the amino acid
    """
    return re.findall(f"[^\W\d]+|\d+|\*", aa)


def split_aa_components(aa):
    """
    Split the PosAA>PosAA string into a list
    """
    return aa.split('>')


def fix_snv_aa(aa):
    """
    Convert the amino acid change value for `csq` output
    to refPOSalt (e.g. S1188L)

    Arguments:
        aa: amino acide from csq (e.g. 1188S>1188L)
    """
    #aa1, aa2 = split_aa_components(aa)
    if aa:
        aa_comp = split_aa_components(aa)
        aa_ref = list()
        aa_alt = list()
        if len(aa_comp) == 2:
            aa_ref = get_aa_components(aa=aa_comp[0])
            aa_alt = get_aa_components(aa=aa_comp[1])
        else:
            aa_ref = get_aa_components(aa=aa_comp[0])
            aa_alt = get_aa_components(aa=aa_comp[0])
        return ''.join([aa_ref[1], aa_ref[0], aa_alt[1]])
    else:
        return 'NONE'


def fix_del_aa(aa):
    """
    Replace the `csq` provided PosAA>PosAA for deletions
    with the format AAPosAA typical of SNPEff and ANNOVAR
    (e.g. 3674LSGF>3674L to S3675_F3677del)
    """
    var_type = 'del'
    aa1, aa2 = split_aa_components(aa=aa)
    aa_ref = get_aa_components(aa=aa1)
    aa_alt = get_aa_components(aa=aa2)
    aa_alt_length = len(aa_alt[1])
    aa_pos = int(aa_ref[0]) + aa_alt_length
    aa_ref_length = len(aa_ref[1])
    aa_ref_new = re.sub(str(aa_alt[1]), '', aa_ref[1])
    aa_ref_new_pos = int(aa_ref[0]) + int(aa_alt_length)
    aa_ref_new_length = len(aa_ref_new)
    aa_alt_new_pos = aa_ref_new_pos + aa_ref_new_length - 1
    aa_alt_new = str(aa_ref_new)[-1]
    return ''.join([aa_ref_new[0], str(aa_ref_new_pos), '_', aa_alt_new, str(aa_alt_new_pos), var_type])


def fix_ins_aa(aa):
    """

    """
    var_type = 'ins'
    pass


def fix_aa_change(var):
    """
    The `bcftools csq` uses a different form to label
    the amino acid change.  This converts the aa change
    to AAPOSaa

    Arguments:
        var: single VCF variant entry

    Return Values:
        Return a dictionary containing the variant values:
            * chrom
            * pos
            * ref
            * alt
            * Consequence
            * gene
            * protein
            * aa
    """
    aa_change = 'NONE'
    aa = get_aa_from_info(var)
    if var.var_type == 'snp':
        aa_change = fix_snv_aa(aa=aa)
    elif var.var_type == 'indel':
        if var.is_deletion:
            aa_change = fix_del_aa(aa=aa)
        # case of multi-nucleotide polymorphism is still
        # categorized as indel in the VCF entry but should
        # be handled like SNP
        elif len(var.REF) == len(var.ALT[0]):
            aa_change = fix_snv_aa(aa=aa)
        elif re.findall('\*', aa):
            aa_change = aa
        else:
            pass
    return aa_change


def convert_alt_list_to_string(alt):
    """
    The ALT field in the VCF file is represented as a list, can convert this to
    a comma seperated string.
    """
    _vars = list()
    for _var in alt:
        _vars.append(str(_var))
    return ','.join(_vars)


def write_variant_file(vars, out):
    '''
    Write the variant data to a file.
    '''
    with open(out, 'w') as fo_p:
        fo_p.write(_variant_file_header())
        fo_p.write('\n')
        for var in vars:
            if var['gene'] != '' and var['protein'] != '':
                gene_protein = '-'.join([var['gene'], var['protein']])
            else:
                gene_protein = 'NA'
            fo_p.write('\t'.join([var['sample'],
                                  var['chr'],
                                  var['pos'],
                                  var['ref'],
                                  var['alt'],
                                  var['consequence'],
                                  var['gene'],
                                  var['protein'],
                                  gene_protein,
                                  ]))
            fo_p.write('\n')


def _variant_file_header():
    return '\t'.join(['sample',
                      'chr',
                      'pos',
		      'ref',
                      'alt',
                      'Consequence',
                      'gene',
                      'protein',
                      'aa'])


def get_gene_from_bcsq(var):
    """
    Get the gene impacted from the BCSQ string in the VCF
    entry.
    """
    bcsq = str()
    if 'BCSQ' in var.INFO:
        bcsq = var.INFO['BCSQ'][0].split('|')
        return bcsq[1]
    else:
        return 'NONE'


def get_consequence_from_bcsq(var):
    """
    Get the functional consequence from the BCSQ string in the VCF
    entry.
    """
    bcsq = str()
    if 'BCSQ' in var.INFO:
        bcsq = var.INFO['BCSQ'][0].split('|')
        return bcsq[0]
    else:
        return 'NONE'


def get_protein_from_bcsq(var):
    """
    Get the protein impacted from the BCSQ string in the VCF
    entry.
    """
    bcsq = str()
    if 'BCSQ' in var.INFO:
        bcsq = var.INFO['BCSQ'][0].split('|')
        return bcsq[5]
    else:
        return 'NONE'


def process_variant_detail_dictionary(var, sample):
    """
    Construct the variant information for the amino acid table
    in a dictionary
    """
    gene = get_gene_from_bcsq(var=var)
    protein = get_protein_from_bcsq(var=var)
    if not gene:
        gene = ''
    if not protein:
        protein = ''
    var_meta = {'sample' : sample,
                'consequence' : get_consequence_from_bcsq(var=var),
                'gene' : get_gene_from_bcsq(var=var),
                #'protein' : get_protein_from_bcsq(var=var),
                'protein' : fix_aa_change(var=var),
                'chr' : str(var.CHROM),
                'pos' : str(var.POS),
                'ref' : str(var.REF),
                'alt' : str(convert_alt_list_to_string(alt=var.ALT))}
    return var_meta


def main():
    """
    Main program
    """
    args = init_args()
    vars = list()
    vcf_reader = vcf.Reader(filename=args.vcf)
    for var in vcf_reader:
        # only take mutations that have a functional annotation
        if 'BCSQ' not in var.INFO:
            continue
        else:
            _var = process_variant_detail_dictionary(var=var, sample=args.sample)
            vars.append(_var)
    write_variant_file(vars=vars, out=args.output)


if __name__ == '__main__':
    main()


#__END__
