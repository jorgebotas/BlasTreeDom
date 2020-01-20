#!/usr/bin/env python

import argparse
import os
import sys
import time

from datetime import datetime, timedelta
from subprocess import call, PIPE, Popen

from Bio import SeqIO
import pandas as pd

import blast as bl
import file_handler as fh
import genbank_parser as gbp
import graphication as graph
import muscle as ms
import prosite_parser as proparse
import user_interface as ui


def main():
    # e-value threshold for blast
    E_VALUE = '1e-03'
    SEQ_TYPE = 'prot'

    arg_parser = argparse.ArgumentParser(
                 description='BlasTreeDom is a local python package that given \
                 input query FASTA file(s) and genBank file(s) performs: \
                 extraction of CDS information and protein sequences from \
                 input genBank(s), blast analysis (blastp), aligment of \
                 original sequences, computation of \
                 Neighbor-Joining phylogenetic trees using MUSCLE \
                 and mapping of ProSite protein domains, providing a \
                 friendly output through graphication of the results.')
    mutually_exclusive = arg_parser.add_mutually_exclusive_group(required=True)
    mutually_exclusive.add_argument(
                      '-query', type=str, nargs='+',
                      help='FASTA file or directory \
                      containing query protein sequences')
    mutually_exclusive.add_argument(
                      '-ui', action='store_true',
                      help='Launch Friendly User Interfase Mode')
    sseqs = arg_parser.add_mutually_exclusive_group(
                       required=not '-ui' in sys.argv)
    sseqs.add_argument('-genBank', type=str, nargs='+',
                       help='genBank file or directory parsed to extract \
                       CDS protein sequences')
    sseqs.add_argument('-multifasta', type=str, nargs='+',
                       help='FASTA file containing subject protein sequences')
    arg_parser.add_argument('-database', type=str,
                            help='If database already computed, \
                                  it can be provided. However, \
                                  original subject sequences are needed \
                                  (genBank or multifasta)')
    arg_parser.add_argument('-pident', type=float,
               help='Identity percentage threshold for blast analysis')
    arg_parser.add_argument('-cov', type=float,
               help='Coverage (percentage) threshold for blast analysis')
    arg_parser.add_argument('-e_value',type=float,
              help='E-value threshold for blast analysis')
    arg_parser.add_argument('-results_dir', type=str,
              help='Output directory to store results. Default: "results/"')
    arg_parser.add_argument('-graph', action='store_true',
               help='Boolean to graph blast and domains analysis outputs')
    args = arg_parser.parse_args()

    if args.ui:
        args.query, args.genBank, args.multifasta, args.cov, args.pident,\
        args.e_value, args.graph = ui.friendly_user_interfase()

    # Running time starts after getting all input parameters
    time0 = time.time()
    now = str(datetime.now()).rsplit('.', 1)[0]\
                             .replace(' ', '_')\
                             .replace(':', '.')


    if args.results_dir: results = os.path.abspath(args.results_dir)\
                                          .rstrip('/')+'/'
    else: results = 'results/'+now+'/'
    os.makedirs(results, exist_ok=True) # Recursively create results directory

    # Boolean to check whether previous step(s) have been computed
    toBeContinued = False
    logfile = results+'_log'
    open(logfile, 'w').close() # Create logfile
    blast_output = results+'_blast_output'

    # Create a single multifasta and tsv file containing all queries from input
    # (single, several files or directory)
    if len(args.query) == 1:
        query = results+os.path.basename(args.query[0]).rsplit('.', 1)[0]
    else:
        query = results+'query'
    fh.merge_files(
                   input_files=args.query,
                   output_dir=results,
                   output_filename=query+'.fasta'
                   )

    # Create directory for each query, inside results
    fh.fasta2dirs(fasta_file=query+'.fasta', output_dir=results)

    # Create a single multifasta file containing
    # all subject sequences from input (single, several files or directory)
    if not args.genBank:
        if len(args.multifasta) == 1:
            gb_multifasta_filename = results \
                                     + os.path.basename(args.multifasta[0])
        else:
            gb_multifasta_filename = results+'genBank_multifasta.fasta'
        fh.merge_files(
                       input_files=args.multifasta,
                       output_dir=results,
                       output_filename=gb_multifasta_filename
                       )
    else:
        gb_multifasta_filename = results+'genBank_multifasta.fasta'

    if args.database:
        database = args.database
        if not args.genBank and not args.multifasta:
            print('\nOriginal subject sequences through genBank \
                   or FASTA file must be provided\n')
            exit(1)
    else:
        database = results+'database/genBank'

    if args.e_value: e_value = args.e_value
    else: e_value = E_VALUE


    print()
    if args.genBank:
        # Generate combined multifasta with all parsed GenBank files
        print("Generating multifasta from GenBank file(s)")
        gbp.parse_gbs(
                      genBanks=args.genBank,
                      sequence_type=SEQ_TYPE,
                      output_dir=results,
                      output_filename=gb_multifasta_filename
                      )
        toBeContinued = True

    if args.multifasta or toBeContinued:
        # Generate database from created multifasta
        print("Generating database...")
        bl.multifasta2database(
                               multifasta=gb_multifasta_filename,
                               sequence_type=SEQ_TYPE,
                               output_dir=results,
                               output_filename=database,
                               log=logfile
                               )
        toBeContinued = True

    # Perform blastp
    print("Performing blast analysis...")
    bl.blast_compute(
                     query_fasta=query+'.fasta',
                     database_path=database,
                     sequence_type=SEQ_TYPE,
                     e_value=e_value,
                     cov_threshold=args.cov,
                     pident_threshold=args.pident,
                     output_dir=results,
                     output_filename=blast_output,
                     log=logfile
                     )
    # Include complete subject sequences in blast_output
    bl.retrieve_seqs(
                     query_fasta=query+'.fasta',
                     subject_multifasta=gb_multifasta_filename,
                     blast_output=blast_output+'.tsv',
                     output_dir=results,
                     output_filename=blast_output+'.tsv',
                     remove_files=True
                     )

    # Include query_fasta and perform multiple alignment(s) using MUSCLE
    print("Performing multiple alignment(s)...")
    ms.compute_alignments(
                          blast_output=blast_output+'.tsv',
                          output_dir=results
                          )

    # Compute NJ tree using MUSCLE
    print("Computing N-J phylogenetic tree(s)...")
    ms.compute_trees(
                     blast_output=blast_output+'.tsv',
                     output_dir=results,
                     output_filename="NJ.phy",
                     log=logfile
                     )

    # Map domains and store them
    print("Extracting ProSite domains...")
    proparse.find_domains(
                          blast_output=blast_output+'.tsv',
                          output_dir=results,
                          summary=True
                          )

    # Merge blast output, genBank info and ProSite domains (only names)
    # into one tsv file
    fh.merge_results(
                     blast_output=blast_output+'.tsv',
                     genBank_tsv=results+'_genBank_info.tsv',
                     output_dir=results,
                     output_filename='_merged.tsv',
                     include_gb=args.genBank
                     )

    if args.graph:
        print("Creating and storing graphs...")
        graph.blast_plot(
                         blast_output=blast_output+'.tsv',
                         output_dir=results
                         )
        graph.domain_plot(
                          blast_output=blast_output+'.tsv',
                          output_dir=results
                          )

    # Ring bell to notify completion
    print("\nProcess COMPLETED")
    print("Running time: "+ str(timedelta(seconds=(time.time() - time0))))
    call(['echo', '\007'])


    sys.exit()


if __name__ == '__main__':
    main()
