# Privacy Preserving Variant Co-Occurrence Analysis

This repository holds a project by the Privacy Preserving Participation group to attempt to classify variants of unknown significance as either benign or pathogenic, by collecting variant co-occurrence information from a large cohort of participants in a privacy-preserving manner, in the context of the BRCA1 and BRCA2 genes.

## Background Facts

The BRCA1 gene contains 125,951 bases, and produces a protein of 1,863 amino acids. [🔖](http://www.genecards.org/cgi-bin/carddisp.pl?gene=BRCA1) The BRCA2 gene contains 85,405, and produces a protein of 3,418 amino acids. [🔖](http://www.genecards.org/cgi-bin/carddisp.pl?gene=BRCA2) Each individual participant has two distinct copies of the BRCA1 gene and two distinct copies of the BRCA2 gene.

There are 18,188 known variants across the two genes, but only 6,154 have expert classifications. [🔖](http://brcaexchange.org/factsheet)

## Step 1: Prototyping

The first step of the project is to construct a prototype system that generates synthetic participant data, computes co-occurring pairs of variants, applies the [ppp library](https://github.com/ppplab/ppp) to generate privacy-preserving summary statistics about the synthetic participants, and then demonstrates the extent to which useful co-occurrence data can be recovered.

### Data Generation

Synthetic data is generated as follows:

1. A population of a single individual, with a single random base haplotype for each gene, is generated.

2. At each time step, each new individual is produced by creating two haplotypes for each gene. Each created haplotype is formed by recombining the haplotypes present in a randomly chosen parent, and adding novel variants. Each variant is identified by a (position, base) tuple, and is designated as either benign or pathogenic.

3. Individuals carrying two pathogenic variants *in trans* in any gene is removed from the population.

4. A final population of participants is generated. The ground truth for each participant, consisting of the set of variants that they carry on each of their haplotypes, and those variants' classifications, is available.

5. **TODO** A censored observed data set is created, by hiding the true classifications of some portion of variants.

To run the data generation prototype, do something like:

```
[anovak@kolossus cooccur]$ prototype/make-data.py --seed 10 ./test.pop ./test.path
...
Completed generation 997/1000
Completed generation 998/1000
Completed generation 999/1000
Completed generation 1000/1000
Gene 0
        98:C->G: Benign 0.608333333333
        79:C->G: Benign 0.0166666666667
        26:G->A: Pathogenic 0.00833333333333
        4:T->A: Pathogenic 0.00666666666667
        51:C->T: Pathogenic 0.00333333333333
Gene 0: 11/300 individuals carry a pathogenic variant
Successfully saved population
Successfully saved pathogenicity database
```

### Privacy Preserving Analysis

The current prototype implementation simply produces a histogram of co-occurrence pairs before and after the privacy-preserving data transmission, and compares them. The histogram is currently over all cooccurrence tuples that occur; it needs to be over all cooccurrence tuples that are **possible** instead, to reflect the information we will have available.

```
[anovak@kolossus cooccur]$ prototype/process-data.py ./test.pop ./test.path
Computing histograms for Gene 0
Counter({((79, 'G'), (98, 'G')): 9, ((4, 'A'), (98, 'G')): 3, ((26, 'A'), (98, 'G')): 3, ((51, 'T'), (98, 'G')): 2})
[[  3.           4.07093923]
 [  3.          11.60475213]
 [  2.           5.92761797]
 [  9.           7.14245489]]
```

A more complete implementation would include prototype classification logic, based on the following principles:

1. Two pathogenic variants in the same gene will not exist *in trans* (i.e. in opposite gene copies) in the same adult human.

2. If variants A and B co-occur twice, variants A and C co-occur, variants B and D co-occur, (and variants C and D do not co-occur), then it is highly likely that variants A and B have been observed *in trans*. [🔖](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2563222/)



