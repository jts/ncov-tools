#!/usr/bin/env python

"""
A small workaround for building the MN908947.3 database in snpeff.
"""

import os

SNPEFF_PREFIX = 'snpeff-5.0-0'

def main():
    """
    Main function to run shell commands to build the MN908947.3 snpEff database
    """
    snpeff_dir = '/'.join([os.environ['CONDA_PREFIX'], 'share', SNPEFF_PREFIX])
    os.chdir(snpeff_dir)
    os.system('./scripts/buildDbNcbi.sh MN908947.3')


if __name__ == "__main__":
    main()

