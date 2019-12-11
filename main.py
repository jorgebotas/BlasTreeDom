#!/usr/bin/env python

import argparse
import os
import sys

from subprocess import call, PIPE, Popen
from Bio import SeqIO

import genbank_parser as gbp
import blast as bl
import muscle as ms



def main():
    # e-value threshold for blast 
    E_VALUE = "1e-03"

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
        bl.blast_compute(query_fasta=args.query, database_path=database, sequence_type=args.sequence_type, e_value=E_VALUE, output_filename=results+"blast_output", log=logfile, fasta=True)
        # Include query_fasta and perform multiple alignment using MUSCLE
        print("Performing multiple alignment...")
        ms.multiple_alignment(multifasta=results+"blast_output.fasta", query_fasta=args.query, output_filename=alignment, log=logfile)

    if args.alignment or toBeContinued:
        # Compute NJ tree using MUSCLE
        print("Computing NJ phylogenetic tree...")
        ms.compute_NJtree(alignment=alignment, output_filename=results+"NJ.tree", log=logfile)
        print("Process COMPLETED")

    # Ring bell to notify completion
    call(['echo', '\007'])

    sys.exit()


if __name__ == '__main__':
    main()