'''
Suite of tests for the Consensus module
'''

import unittest
import ncov.parser

test_meta = ncov.parser.Meta(file='data/metadata.tsv')
test_meta.import_metadata()

class MetaTest(unittest.TestCase):
    def test_get_meta_for_sample(self):
        sample_meta = test_meta.get_meta_for_sample(sample='sampleA')
        self.assertEqual(sample_meta['qpcr_ct'], '17.4')
        self.assertEqual(sample_meta['collection_date'], '2020-03-02')
        self.assertEqual(sample_meta['num_weeks'], 9)
    def test_get_meta_for_sample_old_date(self):
        sample_meta = test_meta.get_meta_for_sample(sample='sampleC')
        self.assertEqual(sample_meta['collection_date'], '1905-04-01')
        self.assertEqual(sample_meta['num_weeks'], 'NA')
