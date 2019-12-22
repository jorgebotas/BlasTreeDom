#!/usr/bin/env python

import argparse
import os
import sys
import time

import datetime
from subprocess import call, PIPE, Popen
from Bio import SeqIO

import blast as bl
import genbank_parser as gbp
from graphication import blast_grapher, domain_grapher
import muscle as ms
import prosite_parser as proparse



def main():
    # e-value threshold for blast 
    E_VALUE = "1e-03"

    time0 = time.time()

    arg_parser = argparse.ArgumentParser()
    # Create mutually exclusive group of positional arguments
    # Diferent actions will be performed: from only generating tree from alignment to performing the whole script
    mutually_exclusive = arg_parser.add_mutually_exclusive_group(required=True) # One must be provided
    mutually_exclusive.add_argument("-query", type=str)
    mutually_exclusive.add_argument("-alignment", type=str) 
    # Create group for inputs containing query file
    non_aligned = arg_parser.add_mutually_exclusive_group()
    non_aligned.add_argument("-genBank", type=str)
    non_aligned.add_argument("-multifasta", type=str)
    non_aligned.add_argument("-database", type=str)
    # Sequence type argument: either prot or nucl
    arg_parser.add_argument("-sequence_type", choices=['prot', 'nucl'])
    arg_parser.add_argument("-results_dir", type=str)
    arg_parser.add_argument("-graph", action='store_true')
    args = arg_parser.parse_args()


    if args.results_dir: results = args.results_dir
    else: results = "results/"
 
    if args.multifasta: gb_multifasta_filename = args.multifasta
    else: gb_multifasta_filename = results+"genBank_multifasta.fasta"

    if args.database: database = args.database
    else: database = results+"database/genBank"

    if args.alignment: alignment = args.alignment
    else: alignment = results+"alignment.fasta"

    toBeContinued = False # Boolean to check whether previous step(s) have been computed
    logfile = results + 'log'
    blast_output = results + 'blast_output'


    if args.genBank:
        # Generate combined multifasta with all parsed GenBank files
        print("Generating multifasta from GenBank file(s)")
        gbp.gb_dir2multifasta(genBank_dir=args.genBank, sequence_type=args.sequence_type, output_filename=gb_multifasta_filename)
        toBeContinued = True
    
    if args.multifasta or toBeContinued:
        # Generate database from created multifasta
        print("Generating database...")
        bl.multifasta2database(multifasta=gb_multifasta_filename, sequence_type=args.sequence_type, output_filename=database, log=logfile)
        toBeContinued = True

    if args.database or toBeContinued:
        # Perform blastp or blastn for protein or nucleotide sequences respectively
        print("Performing blast analysis...")
        bl.blast_compute(query_fasta=args.query, database_path=database, sequence_type=args.sequence_type, e_value=E_VALUE, output_filename=blast_output, log=logfile, fasta=True)
        # Include query_fasta and perform multiple alignment using MUSCLE
        print("Performing multiple alignment...")
        ms.multiple_alignment(multifasta=blast_output+'.fasta', query_fasta=args.query, output_filename=alignment, log=logfile)

    if args.alignment or toBeContinued:
        # Compute NJ tree using MUSCLE
        print("Computing NJ phylogenetic tree...")
        ms.compute_NJtree(alignment=alignment, output_filename=results+"NJ.tree", log=logfile)

    # Map domains and store them
    print("Extracting ProSite domains...")
    domains_dir = results+'domains/'
    if not os.path.isdir(domains_dir): os.mkdir(domains_dir)
    proparse.extract_domains(input_fasta=blast_output+'.fasta', output_dir=domains_dir)


    if args.graph:
        print("Creating and storing graphs...")
        # Create graphs directory to store plots
        graphs_dir = results+'graphs/'
        if not os.path.isdir(graphs_dir): os.mkdir(graphs_dir)
        blast_grapher.blast_plot(blast_output=blast_output+'.tsv', output_dir=graphs_dir, show=True)
        domain_graphs_dir = "results/graphs/domains/"
        if not os.path.isdir(domain_graphs_dir): os.mkdir(domain_graphs_dir)
        domain_grapher.domain_plot("results/blast_output.tsv", domain_graphs_dir)

    # Ring bell to notify completion
    print("Process COMPLETED")
    print("Total time: {}".format(str(datetime.timedelta(seconds=(time.time() - time0))))) ## CREATE FUNCTION IN SYSTEM MODULE
    call(['echo', '\007'])


    sys.exit()


if __name__ == '__main__':
    main()