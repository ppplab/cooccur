#!/usr/bin/env python
""" 

make-data.py: generate synthetic cooccurrence data for benign and pathogenic
variants

This script will generate synthetic participant data according to an
evolutionary model.

Population data is saved as a pickled list of individuals. Each individual is a
list of genes. Each gene is a list of haplotypes. Each haplotype is a list of
(position, base) variant tuples.

Variant classification data is saved as a pickled list, one per gene, of dicts
from variant tuple to boolean pathogenic flag.

"""
import argparse
import os
import sys
import random

import collections
import itertools
import pickle

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
    
    parser.add_argument("--full_scale", action="store_true",
        help="run at full scale gene/pop size instead of small test size")
    parser.add_argument("--seed", type=int, default=None,
        help="use the given integer seed for the random number generator")    
    parser.add_argument("population_file", type=argparse.FileType("w"),
        help="file to save generated individuals to")     
    parser.add_argument("pathogenicity_file", type=argparse.FileType("w"),
        help="file to save pathogenicity flags to")
        
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

# Here is the code that will generate our synthetic participants

# A participant in the simulator is represented by a list of genes. Each gene
# is a list of two haplotypes. Each haplotype is a sorted list of (position,
# base) tuples, where unspecified positions have the reference base for that
# gene, and specified positions always differ from the reference.

def generate_reference(length):
    """
    Generate a reference string of the given length to be the reference for a gene.
    """
    return "".join((random.choice("ACGT") for i in range(length)))

def recombine(haplotype1, haplotype2, length, p_whole = 0.9):
    """
    Recombine two haplotypes. Each haplotype is a sorted list of (position,
    base) tuples. length gives the length of the gene that the haplotypes
    belong to. p_whole gives the probability of taking entirely one haplotype
    or the other (i.e. the crossover point was outside the gene). Otherwise, a
    crossover point is selected and the haplotypes are crossed (which may or
    may not actually result in a new haplotype).
    """

    if random.random() < p_whole:
        # We should just pick one whole haplotype
        return random.choice([haplotype1, haplotype2])
    else:
        # We need to do a crossover. Pick a point
        crossover_point = random.randint(1, length - 1)

        if random.random() < 0.5:
            # Flip the haplotypes
            temp = haplotype1
            haplotype1 = haplotype2
            haplotype2 = temp

        new_haplotype = ([(pos, base) for pos, base in haplotype1 if pos < crossover_point] +
            [(pos, base) for pos, base in haplotype2 if pos >= crossover_point])

        return new_haplotype

def mutate(haplotype, reference):
    """
    Mutate a haplotype once. Adds a new variant to the haplotype that is not
    already in it. May clobber an existing variant. Allows back-mutation to the
    reference. Returns the modified haplotype.
    """
    
    # We need a dict from position to mutation, if any
    hap_dict = dict(haplotype)

    # Where will we mutate
    mutation_pos = random.randint(0, len(reference) - 1)
    # What does the sample have there?
    sample_base = hap_dict.get(mutation_pos, reference[mutation_pos])
    # What will we adopt instead?
    mutation_base = random.choice([base for base in ["A", "C", "G", "T"] if base != sample_base])

    # Construct a mutation record
    mutation = (mutation_pos, mutation_base) if mutation_base != reference[mutation_pos] else None

    # Take the part of the haplotype before the mutation
    new_haplotype = [(pos, base) for pos, base in haplotype if pos < mutation_pos]

    if mutation is not None:
        # We have a non-reference mutation to inject
        new_haplotype.append(mutation)

    # Now add the trailing part of the haplotype
    new_haplotype += [(pos, base) for pos, base in haplotype if pos > mutation_pos]

    return new_haplotype

def decide_if_pathogenic(variant, database, p_pathogenic = 0.5):
    """
    Decide if the given (pos, base) variant is pathogenic in the given database
    dict. If it is not in the database, decides whether it will be pathogenic
    or not and adds it to the database. Each gene must have its own database.
    """

    # The database holds True if pathogenic and False otherwise

    if not database.has_key(variant):
        # We need to decide on an effect
        database[variant] = random.random() < p_pathogenic
    return database[variant]

def has_pathogenic_in_trans(haplotype_pair, database):
    """
    Given a pair of haplotypes for a gene and that gene's pathogenicity
    database, return True if there are pathogenic variants in trans, and false
    otherwise.
    """
    
    # We will set a flag if we have a pathogenic variant in heach haplotype
    pathogenic1 = False
    pathogenic2 = False
    
    for variant in haplotype_pair[0]:
        if decide_if_pathogenic(variant, database):
            pathogenic1 = True
            break
    for variant in haplotype_pair[1]:
        if decide_if_pathogenic(variant, database):
            pathogenic2 = True
            break
            
    return pathogenic1 and pathogenic2

def make_next_generation(population, references, databases, pop_size, p_mutation = 0.001):
    """
    Make the next generation population from the previous generation
    population.  Also takes a list of reference strings for all genes, and a
    list of pathogenicity databases for all genes.  The population will be the
    given target size before accounting for pathogenicity in offspring.
    
    p_mutation gives the probability of having another mutation in a copy of a
    gene; mutations are exponentially distributed.
    """

    new_population = []

    for i in range(pop_size):
        # Make each child
        
        # Pick two parents (may be the same person)
        parent1 = random.choice(population)
        parent2 = random.choice(population)
        
        # Define the offspring (list of genes, which are lists of haplotypes)
        offspring = []

        for gene in range(len(references)):
            # For each gene

            # Decide what haplotypes are sent from the parents
            contributed1 = recombine(parent1[gene][0], parent1[gene][1], len(references[gene]))
            contributed2 = recombine(parent2[gene][0], parent2[gene][1], len(references[gene]))
            
            # Do mutations
            mutation_count = 0
            while random.random() < p_mutation:
                # Mutations are exponential
                contributed1 = mutate(contributed1, references[gene])
                mutation_count += 1
            while random.random() < p_mutation:
                # Do the other copy
                mutation_count += 1
                contributed2 = mutate(contributed2, references[gene])
                
            # We won't decide on pathogenicity until later, so we don't care what mutation is made.
            
            # Append gene to offspring
            offspring.append([contributed1, contributed2])
        
        # Put the finished offspring in the population
        new_population.append(offspring)
        
    # Now cull the population by pathogenicity.  Maybe we should have done this
    # earlier but it sort of makes more sense to do it after all the offspring
    # are made.

    # Who has no pathogenic variants in trans?
    survivors = []
    
    for genome in new_population:
        # For each diploid offspring generated
        
        # We will set this flag if we find pathogenic variants in trans
        pathogenic_in_trans = False
        
        for gene in range(len(references)):
            # For each gene in the simulation
            if has_pathogenic_in_trans(genome[gene], databases[gene]):
                # If this diploid genome has pathogenic variants in trans, it
                # will not make it to the next generation
                pathogenic_in_trans = True
                break
                
        if not pathogenic_in_trans:
            survivors.append(genome)
            
    # TODO: Note that it is technically possible for everyone to die!
    return survivors
    
def compute_frequencies(population):
    """
    Given a population, consisting of a list of individuals, who are lists of
    genes, which are haplotype lists of variants, compute variant frequencies.
    
    Returns a list of dicts from variant to frequency in the haplotypes in the
    population.
    
    Population must not be empty.
    """
    
    # How many genes are there
    num_genes = len(population[0])
    
    # This will be filled with the frequency dicts
    frequency_dicts = []
    
    for gene in range(num_genes):
        # For each gene
        variant_counts = collections.Counter()
        
        for individual in population:
            # For each individual
            for haplotype in individual[gene]:
                # For each haplotype they have of this gene
                for variant in haplotype:
                    # For each variant, count it
                    variant_counts[variant] += 1
            
        # Compute the frequency dict
        frequencies = {variant: float(count)/(len(population) * 2) for variant, count in variant_counts.iteritems()}
        
        # And save it
        frequency_dicts.append(frequencies)
    
    return frequency_dicts
            
    
def evolve_population(gene_lengths, pop_size, generations):
    """
    Evolve a population with genes of the specified lengths, over the specified
    number of generations, with a target population size of the specified
    value.
    
    Returns the final population, the gene references, and the pathogenicity
    database list per gene.
    """
    
    # Make the references
    references = [generate_reference(length) for length in gene_lengths]
    
    # And the pathogenicity databases per gene
    databases = [dict() for length in gene_lengths]
    
    # Make a single individual
    founder = []
    for length in gene_lengths:
        # Put an empty pair of haplotypes for each gene
        founder.append([[], []])
    
    # Make that the only person in the population
    population = [founder]
    
    for i in range(generations):
        # Run a bunch of generations
        population = make_next_generation(population, references, databases, pop_size)
        
        print("Completed generation {}/{}".format(i + 1, generations))
        
    # Return the final population
    return (population, references, databases)
    
    
def report_frequencies(population, references, databases, min_threshold = 0.0001):
    """
    Print out a report on allele frequency and pathogenicity status
    """
    
    # Compute allele frequencies
    frequency_dicts = compute_frequencies(population)
    
    for gene, frequencies in enumerate(frequency_dicts):
        # For each gene
        print("Gene {}".format(gene))
        for variant, frequency in sorted(frequencies.iteritems(), reverse=True, key=lambda (first, second): second):
            # For each variant, from common to rare
            if frequency < min_threshold:
                # Unless it's too rare
                return
            # Print its frequency and status
            is_pathogenic = decide_if_pathogenic(variant, databases[gene])
            print("\t{}:{}->{}: {} {}".format(variant[0], references[gene][variant[0]], variant[1],
                "Pathogenic" if is_pathogenic else "Benign", frequency))
                
def report_genotypes(population, databases):
    """
    Print out a report on what fraction of individuals have any pathogenic variant in each gene.
    """
    
    for gene, database in enumerate(databases):
        # For each gene
        
        # Who has any pathogenic alleles?
        pathogenic_count = 0
        
        for individual in population:
            # Does this individual have a pathoigenic allele?
            has_pathogenic = False
            for variant in individual[gene][0]:
                # For each variant in allele 0
                if decide_if_pathogenic(variant, database):
                    # Note if it's pathogenic
                    has_pathogenic = True
                    break
            for variant in individual[gene][1]:
                # For each variant in allele 1
                if decide_if_pathogenic(variant, database):
                    # Note if it's pathogenic
                    has_pathogenic = True
                    break
            
            if has_pathogenic:
                # Count up individuals with a pathogenic variant
                pathogenic_count += 1
            
        print("Gene {}: {}/{} individuals carry a pathogenic variant".format(gene, pathogenic_count, len(population)))
                
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    if options.seed is not None:
        # Seed the RNG for repeatability
        random.seed(options.seed)
    
    if options.full_scale:
        # Run the simulation at real scale
        # TODO: this takes an hour and we really should cache the simulated data
        population, references, databases = evolve_population([125951, 85405], 30000, 1000)
    else:
        # Run the simulation at a smaller scale
        population, references, databases = evolve_population([100], 300, 1000)
    
    
    # Report the final variant frequencies
    report_frequencies(population, references, databases)
    
    # And the overall summary statistics
    report_genotypes(population, databases)
    
    # Serialize the simulated data as pickled individuals. Each individual is a
    # list of genes, with two haplotypes for each gene, and a list of variant
    # tuples as each haplotype.
    pickle.dump(population, options.population_file)
    print("Successfully saved population")
    
    # Serialize the pathogenicity dicts, one per gene, from variant tuple to
    # pathogenicity determination
    pickle.dump(databases, options.pathogenicity_file)
    print("Successfully saved pathogenicity database")
    
    
    
    
    
if __name__ == "__main__" :
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
        

