'''
A module for handling the metadata file.
'''

import os
import sys
import re
import csv
import datetime

class Meta():
    '''
    The Meta class for handling the metadata.tsv file.
    '''

    def __init__(self, file, start_date='2020-01-01', delimiter='\t'):
        '''
        Initialize the Meta object.
        '''
        if os.path.exists(file):
            self.file = file
        else:
            return None
        self.start_date = start_date
        self.delimiter = delimiter
        self.data = dict()


    def import_metadata(self, sample_id='sample', ct_id='ct', date_id='date'):
        '''
        Import the metadata file and obtain the sample, ct and collection date.

        Arguments:
            * sample_id:    the column label representing the sample name
                            (default: 'sample')
            * ct_id:        the column label representing the ct value
                            (default: 'ct')
            * date:         the column label representing the date value
                            (default: 'date')
            * delimiter:    the file delimiter in the text file (default: '\t')

        Return Value:
            Returns a dictionary containing metadata for all samples
        '''
        try:
            with open(self.file) as file_p:
                data = dict()
                meta_reader = csv.DictReader(file_p, delimiter=self.delimiter)
                for line in meta_reader:
                    if line[date_id] != 'NA':
                        num_months = get_number_of_months(date=line[date_id],
                                                          start_date=self.start_date)
                        num_weeks = get_number_of_weeks(date=line[date_id],
                                                        start_date=self.start_date)
                        if num_weeks <= 0:
                            num_weeks = 'NA'
                    else:
                        num_months = 'NA'
                        num_weeks = 'NA'
                    data[line[sample_id]] = {'qpcr_ct': line[ct_id],
                                             'collection_date': line[date_id],
                                             'num_months' : num_months,
                                             'num_weeks' : num_weeks}
                self.data = data
                return data
        except:
            print("Invalid metadata file")


    def get_meta_for_sample(self, sample):
        '''
        Get the meta information for a single sample
        '''
        if sample in self.data:
            return self.data[sample]


def get_number_of_months(date, start_date='2020-01-01', date_format='%Y-%m-%d'):
    '''
    Get the number of months from the collection date to a start date.
    '''
    required_date_format = '%Y-%m-%d'
    end_date = datetime.datetime.strptime(date, required_date_format)
    start_date = datetime.datetime.strptime(start_date, date_format)
    num_months = (end_date.year - start_date.year) * 12 + (end_date.month - start_date.month)
    return num_months


def get_number_of_weeks(date, start_date='2020-01-01', date_format='%Y-%m-%d'):
    '''
    Get the number of weeks from the collection date to a start date.
    '''
    required_date_format = '%Y-%m-%d'
    end_date = datetime.datetime.strptime(date, required_date_format)
    st_date = datetime.datetime.strptime(start_date, date_format)
    num_weeks = (end_date.year - st_date.year) * 52
    end_week = int(datetime.datetime.strptime(date, required_date_format).strftime('%V'))
    start_week = int(datetime.datetime.strptime(start_date, date_format).strftime('%V'))
    num_weeks = num_weeks + (end_week - start_week)
    return num_weeks
