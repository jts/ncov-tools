'''
Functions to handle and process primers.
'''

import os
import sys
import re
import pandas as pd
import csv

def import_bed_file(bed):
    '''
    Import a BED file containing regions and return a matrix containing those
    regions.

    Arguments:
        * bed: full path to a BED file to process

    Return Value:
        Returns a list of primers from the BED file
    '''
    primers = []
    with open(bed, 'r') as bed_f:
        for line in bed_f:
            line = line.rstrip()
            tmp_data = line.split(sep='\t')
            primers.append(tmp_data)
    return primers


def create_primer_pairs(primers, left='_LEFT.*', right='_RIGHT.*'):
    '''
    From a list of primers, create a dictionary containing the left and right
    regions in a single entry.

    Arguments:
        * primers:  a list containing primers to process
        * left:     string pattern for the left primer in the id
        * right:    string pattern for the right primer in the id

    Return Value:
        Returns a dictionary with the primer ID as key and the following:
            * left_start:       start position of the left primer
            * left_end:         end position of the left primer
            * right_start:      start position of the right primer
            * right_end:        end position of the right primer
    '''
    primer_pairs = {}
    pattern = ''.join(["(", left, '|', right, ")"])

    for primer in primers:
        primer_pool = primer['PoolName']
        primer_id = re.sub(pattern, '', primer['Primer_ID'])
        if primer_id not in primer_pairs.keys():
            primer_pairs[primer_id] = dict()

        if re.search(left, primer['Primer_ID']):
            primer_pairs[primer_id]['left_name'] = primer['Primer_ID']
            primer_pairs[primer_id]['left_start'] = int(primer['start'])
            primer_pairs[primer_id]['left_end'] = int(primer['end'])
            primer_pairs[primer_id]['primer_pool'] = primer['PoolName']
        elif re.search(right, primer['Primer_ID']):
            primer_pairs[primer_id]['right_name'] = primer['Primer_ID']
            primer_pairs[primer_id]['right_start'] = int(primer['start'])
            primer_pairs[primer_id]['right_end'] = int(primer['end'])
            primer_pairs[primer_id]['primer_pool'] = primer['PoolName']
        else:
            print('skipping')
            print(primer)
    return primer_pairs

def create_amplicon_range(primer_pairs, pattern):
    '''
    Convert the primer pairs to a BED formatted amplicon range list.

    Arguments:
        * primer_pairs:     output from the create_primer_pairs function
        * pattern:          pattern in the amplicon key to get the amplicon
                            number

    Return Values:
        Returns a list of amplicon ranges in BED format
    '''
    amplicon_range = []
    for amplicon in primer_pairs:
        if all (k in primer_pairs[amplicon] for k in ('left_start', 'right_end')):
            amplicon_num = re.sub(pattern, '', amplicon)
            amplicon_range.append(create_bed_item(ref='MN908947.3',
                                                start=primer_pairs[amplicon]['left_end'],
                                                end=primer_pairs[amplicon]['right_start'],
                                                name=amplicon_num,
                                                score=primer_pairs[amplicon]['primer_pool'],
                                                strand='+'))
    return amplicon_range


def create_bed_item(ref, start, end, name, score, strand):
    '''
    Create a BED formatted list.

    Arguments:
        * ref: reference for genome/chromosome
        * start: start position for amplicon
        * end: end position for amplicon
        * name: name of amplicon
        * score: score of amplicon
        * strand: strand of the amplicon

    Return Values:
        Returns a list
    '''
    return [str(ref),
            str(start),
            str(end),
            str(name),
            str(score),
            str(strand)]


def create_unique_amplicons(amplicons, offset=30):
    '''
    Create a list of unique regions from the amplicons.

    Arguments:
        * amplicons:    a list of amplicons in BED format from
                        create_amplicon_range
        * offset:       nucleotide count offset (default: 30) 
    '''
    unique_amplicons = []
    for index in range(0, len(amplicons)):
        if index == 0:
            tmp_amplicon = [
                amplicons[index][0],
                amplicons[index][1],
                str(int(amplicons[index+1][1]) - int(offset)),
                amplicons[index][3],
                amplicons[index][4],
                amplicons[index][5]
            ]
        elif index == len(amplicons)-1:
            tmp_amplicon = [
                amplicons[index][0],
                str(int(amplicons[index-1][2]) + int(offset)),
                amplicons[index][2],
                amplicons[index][3],
                amplicons[index][4],
                amplicons[index][5]
            ]
        else:
            tmp_amplicon = [
                amplicons[index][0],
                str(int(amplicons[index-1][2]) + int(offset)),
                str(int(amplicons[index+1][1]) - int(offset)),
                amplicons[index][3],
                amplicons[index][4],
                amplicons[index][5]
            ]
        unique_amplicons.append(tmp_amplicon)
    return unique_amplicons


# the code here to the end of the file was obtained from the ARTIC pipeline
# package:
# https://github.com/artic-network/fieldbioinformatics/blob/master/artic/vcftagprimersites.py

def getPrimerDirection(primerID):
    """Infer the primer direction based on it's ID containing LEFT/RIGHT
    Parameters
    ----------
    primerID : string
        The primer ID from the 4th field of the primer scheme
    """
    if 'LEFT' in primerID:
        return '+'
    elif 'RIGHT':
        return '-'
    else:
        print("LEFT/RIGHT must be specified in Primer ID", file=sys.stderr)
        raise SystemExit(1)


def merge_sites(canonical, alt):
    """Merges a canonical primer site with an alt site, producing an interval that encompasses both
    Parameters
    ----------
    canonical : dict
        The canonical primer site, provided as a dictionary of the bed file row
    alt : dict
        The alt primer site, provided as a dictionary of the bed file row
    Returns
    -------
    dict
        A dictionary of the merged site, where the dict represents a bed file row
    """
    # base the merged site on the canonical
    mergedSite = canonical

    # check the both the canonical and alt are the same direction
    if canonical['direction'] != alt['direction']:
        print(
            "could not merge alt with different orientation to canonical", file=sys.stderr)
        raise SystemExit(1)

    # merge the start/ends of the alt with the canonical to get the largest window possible
    if alt['start'] < canonical['start']:
        mergedSite['start'] = alt['start']
    if alt['end'] > canonical['end']:
        mergedSite['end'] = alt['end']
    return mergedSite


def read_bed_file(fn):
    """Parses a bed file and collapses alts into canonical primer sites
    Parameters
    ----------
    fn : str
        The bedfile to parse
    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - Primer_ID, direction, start, end
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    primers = pd.read_csv(fn, sep='\t', header=None,
                          names=['chrom', 'start', 'end',
                                 'Primer_ID', 'PoolName'],
                          dtype={'chrom': str, 'start': int, 'end': int,
                                 'Primer_ID': str, 'PoolName': str},
                          usecols=(0, 1, 2, 3, 4),
                          skiprows=0)
    if len(primers.index) < 1:
        print("primer scheme file is empty", file=sys.stderr)
        raise SystemExit(1)
    if primers.isnull().sum().sum():
        print("malformed primer scheme file", file=sys.stderr)
        raise SystemExit(1)

    # compute the direction
    primers['direction'] = primers.apply(
        lambda row: getPrimerDirection(row.Primer_ID), axis=1)

    # separate alt primers into a new dataframe
    altFilter = primers['Primer_ID'].str.contains('_alt')
    alts = pd.DataFrame(
        columns=('chrom', 'start', 'end', 'Primer_ID', 'PoolName', 'direction'))
    alts = pd.concat([alts, primers[altFilter]])
    primers = primers.drop(primers[altFilter].index.values)

    # convert the primers dataframe to dictionary, indexed by Primer_ID
    #  - verify_integrity is used to prevent duplicate Primer_IDs being processed
    bedFile = primers.set_index('Primer_ID', drop=False,
                                verify_integrity=True).T.to_dict()

    # if there were no alts, return the bedfile as a list of dicts
    if len(alts.index) == 0:
        return list(bedFile.values())

    # merge alts
    for _, row in alts.iterrows():
        primerID = row['Primer_ID'].split('_alt')[0]

        # check the bedFile if another version of this primer exists
        if primerID not in bedFile:

            # add to the bed file and continue
            bedFile[primerID] = row.to_dict()
            continue

        # otherwise, we've got a primer ID we've already seen so merge the alt
        mergedSite = merge_sites(bedFile[primerID], row)

        # update the bedFile
        bedFile[primerID] = mergedSite

    # return the bedFile as a list
    return [value for value in bedFile.values()]