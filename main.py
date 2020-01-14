#!/usr/bin/env python

import argparse
import os
import sys
import time

from datetime import datetime, timedelta
from subprocess import call, PIPE, Popen
from Bio import SeqIO

import blast as bl
import genbank_parser as gbp
import graphication as graph
import muscle as ms
import prosite_parser as proparse
import file_handler as fh



def main():
    # e-value threshold for blast 
    E_VALUE = '1e-03'
    SEQ_TYPE = 'prot'

    time0 = time.time()
    now = str(datetime.now()).rsplit('.', 1)[0].replace(' ', '_')

    arg_parser = argparse.ArgumentParser()
    # Create mutually exclusive group of positional arguments
    # Diferent actions will be performed: from only generating tree from alignment to performing the whole script
    mutually_exclusive = arg_parser.add_mutually_exclusive_group(required=True) # One must be provided
    mutually_exclusive.add_argument("-query", type=str, help='FASTA file or directory containing query protein sequences') 
    mutually_exclusive.add_argument("-unaligned", type=str, help='FASTA file containing unaligned protein sequences')
    # Create group for inputs containing query file(s)
    non_aligned = arg_parser.add_mutually_exclusive_group()
    non_aligned.add_argument("-genBank", type=str, help='genBank file or directory parsed to extract CDS protein sequences')
    non_aligned.add_argument("-multifasta", type=str, help='FASTA file containing subject protein sequences')
    non_aligned.add_argument("-database", type=str, help='If database already computed, it can be provided. \
                                                        However, original subject sequences are needed (genBank or multifasta)')
    arg_parser.add_argument("-pident", type=float, help='Identity percentage threshold for blast analysis')
    arg_parser.add_argument("-cov", type=float, help='Coverage (percentage) threshold for blast analysis')
    arg_parser.add_argument("-e_value", help='E-value threshold for blast analysis')
    arg_parser.add_argument("-results_dir", type=str, help='Output directory to store results. Default: "results/"')
    arg_parser.add_argument("-graph", action='store_true', help='Boolean to graph blast and domains analysis outputs')
    args = arg_parser.parse_args()


    if args.results_dir: results = args.results_dir
    else: results = 'results/'+now+'/'
 
    if args.multifasta: gb_multifasta_filename = args.multifasta
    else: gb_multifasta_filename = results+'genBank_multifasta.fasta'

    if args.database: database = args.database
    else: database = results+'database/genBank'

    if args.alignment: alignment = args.alignment
    else: alignment = results+'alignment.fasta'

    toBeContinued = False # Boolean to check whether previous step(s) have been computed
    logfile = results+'_log'
    blast_output = results+'_blast_output'

    if not os.path.isdir(results): os.makedirs(results, exist_ok=True)

    # Create a single multifasta and tsv file containing all queries if a directory is provided as argument
    if os.path.isdir(args.query):
        query = results+os.path.basename(args.query)
        fh.merge_files(directory=args.query, output_dir=results, output_filename=query+'.fasta')
    else: query = args.query.rsplit('.', 1)[0]

    # Create directory for each query, inside results
    fh.fasta2dirs(fasta_file=query+'.fasta', output_dir=results)

    if args.genBank:
        # Generate combined multifasta with all parsed GenBank files
        print("Generating multifasta from GenBank file(s)")
        gbp.gb_dir2multifasta(genBank_dir=args.genBank, sequence_type=SEQ_TYPE, output_filename=gb_multifasta_filename)
        toBeContinued = True
    
    if args.multifasta or toBeContinued:
        # Generate database from created multifasta
        print("Generating database...")
        bl.multifasta2database(multifasta=gb_multifasta_filename, sequence_type=SEQ_TYPE, output_dir=results, output_filename=database, log=logfile)
        toBeContinued = True

    if args.database or toBeContinued:
        # Perform blastp or blastn for protein or nucleotide sequences respectively
        print("Performing blast analysis...")
        bl.blast_compute(query_fasta=query+'.fasta', database_path=database, sequence_type=SEQ_TYPE, e_value=E_VALUE,
                         cov_threshold=args.cov, pident_threshold=args.pident, output_filename=blast_output, log=logfile)
        # Include complete subject sequences in blast_output
        bl.retrieve_seqs(query_fasta=query+'.fasta', subject_multifasta=gb_multifasta_filename, blast_output=blast_output+'.tsv', output_dir=results,
                         output_filename=blast_output+'.tsv', remove_files=True)
        toBeContinued = True

    if args.unaligned or toBeContinued:                     
        # Include query_fasta and perform multiple alignment(s) using MUSCLE
        print("Performing multiple alignment(s)...")
        ms.compute_alignments(blast_output=blast_output+'.tsv', output_dir=results)
        toBeContinued = True


    if args.alignment or toBeContinued:
        # Compute NJ tree using MUSCLE
        print("Computing NJ phylogenetic tree(s)...")
        ms.compute_trees(blast_output=blast_output+'.tsv', output_dir=results, output_filename="NJ.tree", log=logfile)

    # Map domains and store them
    print("Extracting ProSite domains...")
    proparse.find_domains(blast_output=blast_output+'.tsv', output_dir=results, summary=True)


    if args.graph:
        print("Creating and storing graphs...")
        graph.blast_plot(blast_output=blast_output+'.tsv', output_dir=results, show=False)
        graph.domain_plot(blast_output=blast_output+'.tsv', output_dir=results, show=False)

    # Ring bell to notify completion
    print("\nProcess COMPLETED")
    print("Running time: "+ str(timedelta(seconds=(time.time() - time0))))
    call(['echo', '\007'])


    sys.exit()


if __name__ == '__main__':
    main()