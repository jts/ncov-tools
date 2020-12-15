#!/usr/bin/env python

'''
Process a directory of SNPEff VCF files and generate a tab-seperated file
containing samples, the functional consequence, the gene and the protein change.
'''

import csv
import sys
import argparse
import re


def parse_args():
    '''
    Process the command line arguments.
    '''
    description = 'Process ANNOVAR data and generate a table of mutations'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f', '--file',
                        help='path to the annotated variant file from ANNOVAR')
    parser.add_argument('-s', '--sample',
                        help='name of the sample to be processed')
    parser.add_argument('-o', '--output',
                        help='filename to write the mutation table')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()



def process_variant_detail_dictionary(var, sample):
    '''
    Process a variant and construct a variant detail dictionary entry. 
    '''
    # var_list = var.INFO['ANN'][0].split('|')
    var_list = var['AAChange.avGene'].split(',')
    var_item = var_list[0].split(':')
    gene = var_item[0]
    protein = ''
    if len(var_item) >= 5:
        protein = re.sub('^p.', '', var_item[4])
    if not gene:
        gene = ''
    if not protein:
        protein = ''
    var_meta = {'sample' : sample,
                'consequence' : var['ExonicFunc.avGene'],
                'gene' : gene,
                'protein' : protein}
    return var_meta


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
                                  var['consequence'],
                                  var['gene'],
                                  var['protein'],
                                  gene_protein]))
            fo_p.write('\n')


def _variant_file_header():
    return '\t'.join(['sample', 'Consequence', 'gene', 'protein', 'aa'])


def main(delimiter='\t'):
    '''
    The main function to execute during application startup.
    '''
    args = parse_args()
    vars = list()
    with open(args.file, 'r') as fp:
        var_reader = csv.DictReader(fp, delimiter=delimiter)
        for var in var_reader:
            _var = process_variant_detail_dictionary(var=var, sample=args.sample)
            vars.append(_var)
    fp.close()
    write_variant_file(vars=vars, out=args.output)


if __name__ == "__main__":
    main()

