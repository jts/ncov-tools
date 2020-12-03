#!/usr/bin/env python

import sys
import csv
import pysam
import argparse
import textwrap as tw


def get_detected_amplicons(file):
    """
    Parse the negative control report and obtain a list of detected amplicons.
    """
    amplicons = set()
    with open(file, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter='\t')
        for record in reader:
            if len(record['amplicons_detected']) > 0:
                _amplicons = record['amplicons_detected'].split(',')
                for amplicon in _amplicons:
                    amplicons.add(amplicon)
    return amplicons


def get_amplicon_dictionary(file, amplicons, column=3, delimiter='_'):
    """
    Create a dictionary record for a set of given amplicons.
    """
    amplicon_dict = dict()
    amplicon_data = list()
    with open(file, 'r') as ifh:
        for amplicon in ifh:
            amplicon = amplicon.strip()
            # amplicon_data.append(amplicon.split('\t'))
            _amplicon_data = amplicon.split('\t')
            _id = get_amplicon_id(amplicon=_amplicon_data[column])
            if _id in amplicons:
                amplicon_dict[_id] = {"start" : _amplicon_data[1], "end" : _amplicon_data[2]}
            else:
                continue
    return amplicon_dict


def get_amplicon_id(amplicon, column=1, delimiter='_'):
    """
    Get the amplicon ID from an amplicon BED entry
    """
    if len(amplicon) > 0:
        amplicon_id = amplicon.split(delimiter)
        return amplicon_id[column]
    else:
        return None


def mask_genome(genome, amplicons, mask='N'):
    """
    Mask a genome FASTA with Ns give a list of positions.
    """
    fasta = pysam.FastxFile(genome)
    for record in fasta:
        sequence = record.sequence
        for _id in amplicons:
            start = int(amplicons[_id]['start']) - 1
            end = int(amplicons[_id]['end']) - 1
            mask_size = end - start + 1
            sequence = sequence[:start] + mask_size * mask + sequence[end+1:]
    return {'header' : record.name, 'sequence' : sequence}


def create_fasta(header, sequence):
    """
    Create a FASTA record (includes header and sequence)
    """
    fasta_record = list()
    fasta_record.append(''.join(['>', header]))
    fasta_record = fasta_record + tw.wrap(str(sequence), width=60)
    return fasta_record


def write_fasta(record, file):
    """
    Write the FASTA sequence to a file
    """
    with open(file, 'w') as ofh:
        for line in record:
            ofh.write(line)
            ofh.write("\n")
    ofh.close()


def main():
    """
    Main method for script
    """
    description = 'Mask amplicons detected in negative controls from a consensus genome'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g', '--genome', help='consensus genome FASTA file to process')
    parser.add_argument('-b', '--bed', help='amplicon BED file')
    parser.add_argument('-n', '--negative_control_report', help='the negative control report')
    parser.add_argument('-o', '--output', help='name of FASTA file to write masked genome to')

    if len(sys.argv) <= 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    detected_amplicons = get_detected_amplicons(file=args.negative_control_report)
    amplicon_dict = get_amplicon_dictionary(file=args.bed, amplicons=detected_amplicons)
    sequence = mask_genome(genome=args.genome, amplicons=amplicon_dict)
    fasta = create_fasta(header=sequence['header'], sequence=sequence['sequence'])
    write_fasta(fasta, args.output)


if __name__ == '__main__':
    main()
