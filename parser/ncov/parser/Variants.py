'''
A Pythone module for handling variant.tsv files from nCoV pipelin.
'''

import os
import sys
import re
import csv

class Variants():
    '''
    The Variants class to handle the variants.tsv file.
    '''

    def __init__(self, file, delimiter='\t'):
        '''
        Initialize the object with the variants.tsv file
        '''
        if os.path.exists(file):
            self.file = file
        else:
            sys.exit("Invalid or missing file.")
        self.delimiter = delimiter


    def get_total_variants(self, indel=True):
        '''
        A method that parses the iVar variants file and returns the total
        number of variants.

        Arguments:
            * indel:    a boolean to determine whether to process indels

        Returns:
            Function returns a dictionary containing the following keys:
                * num_variants_snvs:            total number of variants in the file
                * num_variants_indel:           total number of indels in the file
                * num_variants_indel_triplet:   total number of frameshift indels
        '''
        counter = 0
        counter_snv = 0
        counter_indel = 0
        counter_indel_triplet = 0

        with open(self.file) as file_p:
            variant_reader = csv.DictReader(file_p, delimiter=self.delimiter)
            for data in variant_reader:
                if indel:
                    if len(str(data['ALT'])) > 1:
                        counter += 1
                        counter_indel += 1
                        if is_indel_triplet(data['ALT']):
                            counter_indel_triplet += 1
                    elif len(str(data['ALT'])) == 1:
                        counter += 1
                        counter_snv += 1
                elif not indel:
                    if len(str(data['ALT'])) == 1:
                        counter += 1
                        counter_snv += 1
                    else:
                        continue
        file_p.close()
        return {'num_variants_snvs' : counter_snv,
                'num_variants_indel' : counter_indel,
                'num_variants_indel_triplet' : counter_indel_triplet}

def is_indel_triplet(variant, size=3):
    '''
    Check whether the indel is a multiple of 3bp in size.  We will be using
    this to infer potential
    frameshift indels.

    Arguments:
        * variants: a string represent the variant
        * size:     size of a codon indicating a frameshift

    Return Value:
        Function returns a boolean
    '''
    # remove the leading +/- from the variant
    variant = re.sub('^[+-]', '', variant)
    if len(variant) > 0:
        return not len(variant) % size
