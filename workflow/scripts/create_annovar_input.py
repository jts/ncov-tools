#!/usr/bin/env python
'''
A script to convert variant output from the ncov pipeline.
'''

import os
import os.path
import sys
import csv
import argparse
import vcf
import re


class Variant(object):
    def __init__(self, var):
        '''
        A class to handle variants from the variants.tsv file.
        '''
        self.var = var
    
    def convert_var_to_annovar(self, chr='NC_045512v2'):
        '''
        Convert the var dictionary to ANNOVAR compatible inputs.
        '''
        start = str(int(self.var['POS']) + 1)
        end = str(int(self.var['POS']) + 1)
        if self.is_insertion():
            ref = '-'
            alt = re.sub('^[+]', '', self.var['ALT'])
        elif self.is_deletion():
            ref = re.sub('^[-]', '', self.var['ALT'])
            alt = '-'
            end = int(self.var['POS']) + len(ref)
        else:
            start = str(int(self.var['POS']))
            end = str(int(self.var['POS']))
            ref = self.var['REF']
            alt = self.var['ALT']
            
        return {'chr': chr,
                'start': start,
                'end': end,
                'ref': ref,
                'alt': alt}
    
    def is_insertion(self):
        '''
        check if variant is an insertion
        '''
        return self.var['ALT'].startswith('+')

    def is_deletion(self):
        '''
        check if variant is a deletion
        '''
        return self.var['ALT'].startswith('-')



def convert_vcf_to_annovar(var, chr='NC_045512v2'):
    '''
    Using the vcf object, create the ANNOVAR compatible output as a dictionary.
    '''
    start = 0
    end = 0
    ref = ''
    alt = ''
    if var.is_indel:
        if var.is_deletion:
            ref = re.sub(str(var.ALT[0]), '', var.REF)
            alt = '-'
            start = var.POS + 1
            end = var.POS + len(str(ref))
        else:
            start = var.POS + 1
            end = var.POS + 1
            ref = '-'
            alt = re.sub(var.REF, '', str(var.ALT[0]))
    elif var.is_snp:
        start = var.POS
        end = var.POS
        ref = var.REF
        alt = var.ALT[0]
    else:
        print("invalid variant")
        return
    return {'chr': chr,
            'start': start,
            'end': end,
            'ref': ref,
            'alt': alt}


def main():
    '''
    Main routine for the script
    '''
    parser = argparse.ArgumentParser(description='convert variants to annovar input')
    parser.add_argument('-f', '--file', help='variant file to convert')
    parser.add_argument('-o', '--output', help='output filename')
    parser.add_argument('-d', '--delimiter', default='\t', help='column delimiter for output')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    vars = list()

    # check for filetype
    if args.file.endswith('.variants.tsv'):
        samplename = os.path.basename(args.file).replace('.variants.tsv', '')
        with open(args.file, 'r') as fp:
            var_reader = csv.DictReader(fp, delimiter='\t')
            for var in var_reader:
                ann_var = Variant(var=var)
                vars.append(ann_var.convert_var_to_annovar())
    if args.file.endswith('.pass.vcf') or args.file.endswith('.pass.vcf.gz'):
        if args.file.endswith('.pass.vcf'):
            samplename = os.path.basename(args.file).replace('.pass.vcf', '')
        elif args.file.endswith('.pass.vcf.gz'):
            samplename = os.path.basename(args.file).replace('.pass.vcf.gz', '')
        vcf_reader = vcf.Reader(filename=args.file)
        for var in vcf_reader:
            vars.append(convert_vcf_to_annovar(var=var))

    with open(args.output, 'w') as file_o:
        for var in vars:
            file_o.write(args.delimiter.join([var['chr'],
                                str(var['start']),
                                str(var['end']),
                                str(var['ref']),
                                str(var['alt']),
                                samplename]))
            file_o.write('\n')

if __name__ == '__main__':
    main()
