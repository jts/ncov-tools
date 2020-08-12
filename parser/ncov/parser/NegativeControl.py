'''
A module for processing the NegativeControl output.
'''

import os
import sys
import re
import csv

class NegativeControl(object):
    '''
    A class for handling the output from the negative_control_check.py script.
    '''

    def __init__(self, file, delimiter='\t'):
        '''
        Initialize the object
        '''
        if os.path.exists(file):
            self.file = file
        else:
            sys.exit("Invalid or missing file.")
        self.delimiter = delimiter
    

    def get_control_stats(self, sample):
        '''
        Get a list of stats on the negative controls including:
            * qc - PASS/FAIL
            * genome_covered_bases
            * genome_total_bases
            * genome_coverage_span
            * amplicons_detected
        '''
        data = {'qc': 'NA',
                'genome_covered_bases' : 'NA',
                'genome_total_bases' : 'NA',
                'genome_covered_span' : 'NA',
                'amplicons_detected' : 'NA'}
        with open(self.file, 'r') as file_p:
            csv_reader = csv.DictReader(file_p, delimiter=self.delimiter)
            for record in csv_reader:
                if re.search(sample, record['file']):
                    data[sample] = {'qc' : record['qc'],
                                    'genome_covered_bases' : record['genome_covered_bases'],
                                    'genome_total_bases' : record['genome_total_bases'],
                                    'genome_covered_span' : record['genome_covered_span'],
                                    'amplicons_detected' : record['amplicons_detected']}
        return data
