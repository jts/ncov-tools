#!/usr/bin/env python

import sys
import csv
import pysam
import parasail
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
            _amplicon_data = amplicon.split('\t')
            _id = get_amplicon_id(amplicon=_amplicon_data[column])
            if _id in amplicons:
                # the input BED is 1-based, switch to 0-based here
                amplicon_dict[_id] = {"start" : int(_amplicon_data[1]), "end" : int(_amplicon_data[2])}
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
        len_before = len(sequence)
        for _id in amplicons:
            start = int(amplicons[_id]['start'])
            end = int(amplicons[_id]['end']) - 1
            mask_size = end - start + 1
            sequence = sequence[:start] + mask_size * mask + sequence[end+1:]
        assert(len(sequence) == len_before)
    return {'header' : record.name, 'sequence' : sequence}

def get_sequence(file):
    fasta = pysam.FastxFile(file)
    reference = None
    for record in fasta:
        # For this narrow application there should be exactly one entry in the reference file
        assert(reference is None)
        reference = record
    return reference

def get_alignment(reference_genome, input_genome):
    
    # the dna full matrix supports ambiguity codes, although "N"s are not given free mismatches as we might like
    # the alignments appear good enough for our purpose however
    result = parasail.nw_trace_striped_32(input_genome.sequence, reference_genome.sequence, 10, 1, parasail.dnafull)
    traceback = result.traceback
    columns = 120

    position_map = list()

    reference_index = 0
    input_index = 0

    for (ref, query) in zip(traceback.ref, traceback.query):
        if ref != '-' and query != '-':
            position_map.append( (reference_index, input_index) )
        if ref != '-':
            reference_index += 1
        if query != '-':
            input_index += 1

    return position_map

# translate the coordinates of the amplicon set from reference coordinates
# to the coordinate system of the samples we're masking
def translate_amplicons(position_map, amplicon_dict):

    out_amplicons = dict()

    for amplicon_id, amplicon in amplicon_dict.items():

        # Search the position map for the min/max base in the input genome
        # that is mapped to a base within this amplicon
        reference_start = amplicon['start']

        # we set the reference end coordinate to be inclusive (within the amplicon)
        # to simplify the logic below and avoid corner cases where there is a deletion
        # around the amplicon boundary
        # when the new amplicon coordinates are set later we set this back to be exclusive
        reference_end = amplicon['end'] - 1

        input_start = None
        input_end = None

        for (reference_position, input_position) in position_map:
            if reference_position >= reference_start and reference_position <= reference_end:

                # this position is within the amplicon of interest
                if input_start is None or input_position < input_start:
                    input_start = input_position
                if input_end is None or input_position > input_end:
                    input_end = input_position
                
        if input_start is not None and input_end is not None:
            out_amplicons[amplicon_id] = { "start":input_start, "end":input_end + 1 }
    return out_amplicons

def create_fasta(header, sequence):
    """
    Create a FASTA record (includes header and sequence)
    """
    fasta_record = list()
    fasta_record.append(''.join(['>', header + "_masked"]))
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
    parser.add_argument('-r', '--reference-genome', help='fasta file containing the reference genome')
    parser.add_argument('-o', '--output', help='name of FASTA file to write masked genome to')
    if len(sys.argv) <= 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    reference_genome = get_sequence(file=args.reference_genome)
    input_genome = get_sequence(file=args.genome)
    
    detected_amplicons = get_detected_amplicons(file=args.negative_control_report)
    amplicon_dict = get_amplicon_dictionary(file=args.bed, amplicons=detected_amplicons)
    
    if len(input_genome.sequence) > 0:
        position_map = get_alignment(reference_genome, input_genome)
        amplicon_dict = translate_amplicons(position_map, amplicon_dict)
    else:
        position_map = list()
        amplicon_dict = dict()

    sequence = mask_genome(genome=args.genome, amplicons=amplicon_dict)
    fasta = create_fasta(header=sequence['header'], sequence=sequence['sequence'])
    write_fasta(fasta, args.output)


if __name__ == '__main__':
    main()
