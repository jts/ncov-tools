#!/usr/bin/env python

"""
A small workaround for building the MN908947.3 database in snpeff.
"""

import os

snpeff_dir_prefix = 'snpeff'

def main():
    """
    Main function to run shell commands to build the MN908947.3 snpEff database
    """
    for dir_name in os.scandir('/'.join([os.environ['CONDA_PREFIX'], 'share'])):
        if dir_name.name.startswith(snpeff_dir_prefix):
            os.chdir('/'.join([os.environ['CONDA_PREFIX'], 'share', dir_name.name]))
            os.system(f'{"/".join([os.environ["CONDA_PREFIX"], "share", dir_name.name, "scripts/buildDbNcbi.sh"])} MN908947.3')


if __name__ == "__main__":
    main()

