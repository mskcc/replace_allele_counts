# replace_allele_counts
This repo helps to replace the allele counts in maf file after genotyping
## Requirements:
- pandas : [v0.16.2](http://pandas.pydata.org/)
- nose : [v1.3.7](http://nose.readthedocs.io/en/latest/)

[![Build Status](https://travis-ci.org/rhshah/replace_allele_counts.svg?branch=master)](https://travis-ci.org/rhshah/replace_allele_counts)
[![codecov](https://codecov.io/gh/rhshah/replace_allele_counts/branch/master/graph/badge.svg)](https://codecov.io/gh/rhshah/replace_allele_counts)

## remove_variants.py
#### Based on Mutation Annotation Format ([MAF](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification))
- Takes in a MAF with FillOut based on GetBaseCountsMultiSample and outputs a MAF (Note: You will loose the comments if any present in the MAF file)

```
python replace_allele_counts.py --help
usage: replace_allele_counts.py [options]

This tool helps to replace the allele counts from the caller with the allele
counts of GetBaseCountMultiSample

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         make lots of noise
  -imaf SomeID.maf, --input-maf SomeID.maf
                        Input maf file which needs to be fixed
  -ifill SomeID.fillout.txt, --fillout SomeID.fillout.txt
                        Input fillout file created by GetBaseCountMultiSample
                        using the input maf
  -omaf SomeID.maf, --output-maf SomeID.maf
                        Output maf file name
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
```