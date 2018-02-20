#!/usr/bin/env python
"""
process-data.py: demonstrate the feasibility of privacy-preserving co-occurrence
analysis for classifying variants

This script will filter the co-occurrence pairs of variants from
a population through a privacy preserving data collection step.

Population data is loaded from a pickled list of individuals. Each individual
is a list of genes. Each gene is a list of haplotypes. Each haplotype is a list
of (position, base) variant tuples.

Variant classification data is loaded from a pickled list, one per gene, of
dicts from variant tuple to boolean pathogenic flag.

"""
import argparse
import os
import sys
import random

import collections
import itertools
import pickle

import numpy
import ppp.prototype

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("population_file", type=argparse.FileType("r"),
        help="file to read individuals from")     
    parser.add_argument("pathogenicity_file", type=argparse.FileType("r"),
        help="file to read pathogenicity flags from")
    parser.add_argument("--epsilon", type=float, default=1.0,
        help="privacy preserving factor")
        
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def make_cooccurrence_lists(population):
    """
    For each gene, output a list, per individual, of all cooccurring pairs of
    variants.
    """
    
    # Make a list for each gene of lists for each person of cooccurring variant pairs
    cooccurring = [[] for gene in population[0]]
    
    for individual in population:
        for gene_num, gene_copies in enumerate(individual):
            # Make a list of cooccurring pairs for this person
            cooccurring_for_person = []
            
            # Compose one big set of all the variants in this person's
            # haplotypes for this gene
            all_variants = set()
            for haplotype in gene_copies:
                for variant in haplotype:
                    # Add each variant the person has in the gene to the set
                    all_variants.add(variant)
                    
            for var1, var2 in itertools.combinations(sorted(all_variants), 2):
                # For each pair of distinct variants in cannonical order
                
                # Record a cooccurrence
                cooccurring_for_person.append((var1, var2))
                
            # Save the person;s results for this gene
            cooccurring[gene_num].append(cooccurring_for_person)
    
    return cooccurring

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Load the population of individuals, and the per-gene variant classifications
    population = pickle.load(options.population_file)
    databases = pickle.load(options.pathogenicity_file)
    
    # Construct the big cooccurrence lists of all cooccurrence instances, by gene
    cooccurring = make_cooccurrence_lists(population)
    
    # Make a PPP thing to privacy-preserve the data
    preserver = ppp.prototype.Prototype(options.epsilon)
    
    for gene_num, gene_cooccurrences in enumerate(cooccurring):
        # For each gene and the data for all the people for that gene
        print("Computing histograms for Gene {}".format(gene_num))
        
        # Compute a histogram with a counter to make sure we understand the PPP histograms
        cooccur_counts = collections.Counter()
        for person in gene_cooccurrences:
            for cooccurrence in person:
                cooccur_counts[cooccurrence] += 1
        print cooccur_counts
        
        # Flatten cooccurrences all into one list.
        # Also reduce tuples to strings because PPP thinks the data is the wrong shape and the tuples aren't symbols
        # TODO: Make PPP understand multiple data items per participant and tuples as data items
        cooccurrences_flat = [str(co) for person in gene_cooccurrences for co in person]
        
        # Compute true and private histograms
        true_histogram = preserver.get_histogram(cooccurrences_flat, False)
        private_histogram = preserver.get_histogram(cooccurrences_flat, True)
        
        # Make into 2d arrays
        # TODO: Require numpy 1.10 and stack
        true_histogram = numpy.reshape(true_histogram, (len(true_histogram), 1))
        private_histogram = numpy.reshape(private_histogram, (len(private_histogram), 1)) 
        
        # Compare side by side
        both_histograms = numpy.concatenate((true_histogram, private_histogram), axis=1)
        print(both_histograms)
    
if __name__ == "__main__" :
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
        

