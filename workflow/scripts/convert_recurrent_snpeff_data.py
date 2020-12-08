#!/usr/bin/env python

'''
Process a directory of SNPEff VCF files and generate a tab-seperated file
containing samples, the functional consequence, the gene and the protein change.
'''


import vcf
import sys
import argparse
import re

sequence_ontology_term = [
    'chromosome_number_variation',
    'exon_loss_variant',
    'frameshift_variant',
    'rare_amino_acid_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'start_lost',
    'stop_gained',
    'transcript_ablation',
    '3_prime_UTR_truncation',
    'exon_loss',
    '5_prime_UTR_truncation',
    'exon_loss_variant',
    'coding_sequence_variant',
    'conservative_inframe_deletion',
    'conservative_inframe_insertion',
    'disruptive_inframe_deletion',
    'disruptive_inframe_insertion',
    'missense_variant',
    'regulatory_region_ablation',
    'splice_region_variant',
    'TFBS_ablation',
    '5_prime_UTR_premature_start_codon_gain_variant',
    'initiator_codon_variant',
    'splice_region_variant',
    'start_retained',
    'stop_retained_variant',
    'synonymous_variant',
    '3_prime_UTR_variant',
    '5_prime_UTR_variant',
    'coding_sequence_variant',
    'conserved_intergenic_variant',
    'conserved_intron_variant',
    'downstream_gene_variant',
    'exon_variant',
    'feature_elongation',
    'feature_truncation',
    'gene_variant',
    'intergenic_region',
    'intragenic_variant',
    'intron_variant',
    'mature_miRNA_variant',
    'miRNA',
    'NMD_transcript_variant',
    'non_coding_transcript_exon_variant',
    'non_coding_transcript_variant',
    'regulatory_region_amplification',
    'regulatory_region_variant',
    'TF_binding_site_variant',
    'TFBS_amplification',
    'transcript_amplification',
    'transcript_variant',
    'upstream_gene_variant']


def parse_args():
    '''
    Process the command line arguments.
    '''
    description = 'Process SNPEff data and generate a table of recurrent mutations'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f', '--file',
                        help='path to the annotated VCF file from SNPEff')
    parser.add_argument('-s', '--sample',
                        help='name of the sample to be processed')
    parser.add_argument('-o', '--output',
                        help='filename to write the mutation table')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()



def create_gene_mutation(gene, variant):
    '''
    Create a string with gene-protein_change 
    '''
    return '-'.join([gene, variant])


def process_variant_detail_dictionary(var, sample):
    '''
    Process a variant and construct a variant detail dictionary entry. 
    '''
    var_list = var.INFO['ANN'][0].split('|')
    gene = var_list[3]
    protein = re.sub('^p.', '', var_list[10])
    if not gene:
        gene = ''
    if not protein:
        protein = ''
    var_meta = {'sample' : sample,
                'consequence' : var_list[1],
                'gene' : gene,
                'protein' : protein,
                'chr' : str(var.CHROM),
                'pos' : str(var.POS),
                'ref' : str(var.REF),
                'alt' : str(var.ALT)}
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
                                  gene_protein,
                                  var['chr'],
                                  var['pos'],
                                  var['ref'],
                                  var['alt']]))
            fo_p.write('\n')


def _variant_file_header():
    return '\t'.join(['sample', 'Consequence', 'gene', 'protein', 'aa', 'chr', 'pos', 'ref', 'alt'])


def main():
    '''
    The main function to execute during application startup.
    '''
    args = parse_args()
    vcf_reader = vcf.Reader(filename=args.file)
    vars = list()
    for var in vcf_reader:
        _var = process_variant_detail_dictionary(var=var, sample=args.sample)
        vars.append(_var)
    write_variant_file(vars=vars, out=args.output)


if __name__ == "__main__":
    main()

