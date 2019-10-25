# Genomic Interval Similarity

Given a set of BED files, calculate the similarity between the intervals, and optionally perform _a priori_ or _de novo_ differential analysis.

## Usage

```shell
gis [-h] [-n NAMES] [-o PREFIX] bed [bed ...]

positional arguments:
  bed                   BED file(s)

optional arguments:
  -h, --help            show this help message and exit
  -n NAMES, --names NAMES
                        Comma-separated labels for BED files (default: None)
  -o PREFIX, --prefix PREFIX
                        Prefix for output files. (default: similarity)
```

## Installation

```shell
git clone https://github.com/jrhawley/genomic-interval-similarity.git
cd genomic-interval-similarity/
python setup.py install
```
