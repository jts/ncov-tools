'''
A Python module for handling variants from a VCF file from the ONT platform.
'''

import os
import sys
import vcf

class Vcf():
    '''
    A class for handling VCF files from the ONT platform.
    '''

    def __init__(self, file):
        '''
        Initialize the CovVcf object.
        '''
        if os.path.exists(file):
            self.file = file
        else:
            sys.exit("Invalid or missing file.")


    def get_variant_counts(self):
        '''
        Get the count of snvs and indels from the vcf file.
        '''
        counter = 0
        counter_snv = 0
        counter_indel = 0
        counter_indel_triplet = 0
        var_dict = dict()
        vcf_reader = vcf.Reader(filename=self.file)
        for var in vcf_reader:
            var_id = create_vcf_variant_id(var=var)
            if var_id not in var_dict:
                var_dict[var_id] = 1
            counter += 1
            if var.is_indel:
                counter_indel += 1
                if is_indel_triplet(ref=var.REF, alt=var.ALT[0]):
                    counter_indel_triplet += 1
                else:
                    continue
            if var.is_snp and (len(var.ALT[0]) == 1):
                counter_snv += 1
        return {'num_variants_snvs' : counter_snv,
                'num_variants_indel' : counter_indel,
                'num_variants_indel_triplet' : counter_indel_triplet}


def is_indel_triplet(ref, alt, size=3):
    '''
    Check whether the indel is 3bp in size.  We will be using this to infer potential
    frameshift indels.

    Arguments:
        * ref:      a string representing the reference alleles
        * alt:      a string representing the alternate alleles
        * size:     size of a codon indicating a frameshift

    Return Value:
        Function returns a boolean
    '''
    # remove the leading +/- from the variant
    indel_length = abs(len(str(ref)) - len(str(alt)))
    if indel_length > 0:
        return not indel_length % size


def create_vcf_variant_id(var):
    '''
    Create a variant ID consisting of position, reference allele,
    alternate allele as a concatenated string.  Note that the alternate
    allele is provided as a list and needs to be set as a string
    and merged.

    Arguments:
        * var: a VCF record object from pyvcf

    Return Value:
        Returns a string containing a concatenation of variant
        object identifiers.
    '''
    var_alt = []
    for _var in var.ALT:
        var_alt.append(str(_var))
    var_alt_str = ''.join(var_alt)
    return '_'.join([str(var.POS), var.REF, var_alt_str])
