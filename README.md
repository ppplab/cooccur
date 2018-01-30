# Privacy Preserving Variant Co-Occurrence Analysis

This repository holds a project by the Privacy Preserving Participation group to attempt to classify variants of unknown significance as either benign or pathogenic, by collecting variant co-occurrence information from a large cohort of participants in a privacy-preserving manner, in the context of the BRCA1 and BRCA2 genes.

## Background Facts

The BRCA1 gene contains 125,951 bases, and produces a protein of 1,863 amino acids. [🔖](http://www.genecards.org/cgi-bin/carddisp.pl?gene=BRCA1) The BRCA2 gene contains 85,405, and produces a protein of 3,418 amino acids. [🔖](http://www.genecards.org/cgi-bin/carddisp.pl?gene=BRCA2) Each individual participant has two distinct copies of the BRCA1 gene and two distinct copies of the BRCA2 gene.

There are 18,188 known variants across the two genes, but only 6,154 have expert classifications. [🔖](http://brcaexchange.org/factsheet)

## Step 1: Prototyping

The first step of the project will be to construct a prototype system that will generate synthetic participant data, compute co-occurring pairs of variants, apply the [ppp library](https://github.com/ppplab/ppp) to generate privacy-preserving summary statistics about the synthetic participants, and then demonstrate the extent to which useful co-occurrence data can be recovered.

The most trivial implementation of this prototype will simply produce a histogram of co-occurrence pairs before and after the privacy-preserving data transmission, and compare them. A more complete implementation would include prototype classification logic, based on the following principles:

1. Two pathogenic variants in the same gene will not exist *in trans* (i.e. in opposite gene copies) in the same adult human.

2. If variants A and B co-occur twice, variants A and C co-occur, variants B and D co-occur, (and variants C and D do not co-occur), then it is highly likely that variants A and B have been observed *in trans*. [🔖](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2563222/)


