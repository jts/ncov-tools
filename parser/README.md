# ncov_parser

[![Build Status](https://travis-ci.com/rdeborja/ncov_parser.svg?branch=master)](https://travis-ci.com/rdeborja/ncov_parser) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The `ncov_parser` package provides a suite of tools to parse the files generated
in the Nextflow workflow and provide a QC summary file.  The package requires
several files including:
* <sample>.variants.tsv/<sample>.pass.vcf
* <sample>.per_base_coverage.bed
* <sample>.primertrimmed.consensus.fa
* alleles.tsv

An optional metadata file with qPCR ct and collection date values can be
included.

In addition, `bedtools` should be run to generate a
`<sample>.per_base_coverage.bed` file to generate mean and median depth of
coverage statistics.


## Installation
After downloading the repository, the package can be installed using `pip`:
```
git clone git@github.com:jts/ncov-tools.git
cd ncov-tools/parser
pip install .
```


## Usage
The library consists of several functions that can be imported.
```
import ncov.parser
```
Several classes are available representing the different files that can
be processed.
```
ncov.parser.Alleles
ncov.parser.Consensus
ncov.parser.Meta
ncov.parser.PerBaseCoverage
ncov.parser.Variants
ncov.parser.Vcf
```

Similarly, wrapper scripts for creating a standard format output can be found in
`ncov.parser.qc`
```
import ncov.parser.qc as qc
qc.write_qc_summary_header()
qc.write_qc_summary()
```

### Top levels scripts
In the `bin` directory, several wrapper scripts exist to assist in generating
QC metrics.

To create sample level summary qc files, use the `get_qc.py` script:
```
get_qc_summary.py --variants <sample>.variants.tsv or <sample>.pass.vcf
--coverage <sample>.per_base_coverage.bed --meta <metadata>.tsv
--consensus <sample>.primertrimmed.consensus.fa [--indel]
```

Note the `--indel` flag should only be present if indels will be used in the
calculation of variants.

Once this is complete, we can use the `collect_qc_summary.py` script to
aggregate the sample level summary files into a single run tab-separate file.
```
collect_qc_summary.py --path <path to sample.summary.qc.tsv files>
```

Note that this tool has been used in conjunction with the [@jts `ncov-tools`](https://github.com/jts/ncov-tools)
suite of tools. 

## License
MIT
