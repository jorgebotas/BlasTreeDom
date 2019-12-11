#!/usr/bin/env python

import argparse
import os
import sys

import numpy as np
import pandas as pd
from subprocess import call, PIPE, Popen


def multifasta2database(multifasta, sequence_type, output_filename="database/subject", log='/dev/null'):
    """ Create database from multifasta using makeblastdb """
    # Create directory for database files if it does not already exist
    basename = os.path.basename(output_filename)
    database_path = output_filename.replace(basename, "")
    try:
        os.mkdir(database_path)
    except:
        pass
    call(['makeblastdb', '-in', multifasta, '-dbtype', sequence_type, '-logfile' , log, '-out', output_filename])
    return 


def save_multifasta(input_file = "blast_output.tsv", output_filename = "blast_output.fasta"):
    blast_output = pd.read_csv(input_file, delimiter='\t')
    output_file = open(output_filename, 'w')
    for dummy_idx, hit in blast_output.iterrows():
        output_file.write(">" + str(hit.sseqid) + "\n")
        output_file.write(str(hit.sseq) + "\n")
    output_file.close
    return



def blast_compute(query_fasta, database_path, sequence_type, e_value,  output_filename = "blast_output", log='/dev/null', fasta = False, headers = True):
    """ Perform blastp or blastn analysis for protein or nucleotide sequences respectively """
    if sequence_type == "prot":
        blast_type = 'blastp'
    elif sequence_type == "nucl":
        blast_type = 'blastn'

    outfmt = '6 qseqid sseqid qcovs qstart qend pident evalue sseq'
    call([blast_type, '-query', query_fasta, '-db', database_path, '-evalue', e_value, '-out', output_filename + '.tsv', '-logfile', log, '-outfmt', outfmt])

    if headers:
        # Include headers in output .tsv file
        first_line = outfmt[2:].replace(" ", "\t")
        blast = open(output_filename + '.tsv', 'r')
        blast_output = blast.read()
        blast.close()
        output = open(output_filename + '.tsv', 'w')
        output.write(first_line + '\n')
        output.write(blast_output)
        output.close()
    if fasta:
        save_multifasta(output_filename + '.tsv', output_filename + '.fasta')
    return


def main():
    # e-value threshold 
    E_VALUE = "1e-03" 

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-query', type=str, required=True)
    exclusive = arg_parser.add_mutually_exclusive_group(required=True)
    exclusive.add_argument('-subject', type=str)
    exclusive.add_argument('-database', type=str)
    arg_parser.add_argument("-sequence_type", choices=['prot', 'nucl'], required=True)
    arg_parser.add_argument('-fasta', action='store_true')
    arg_parser.add_argument('-headers', action='store_true')
    args = arg_parser.parse_args()

    if args.database: 
        database = args.database
    else: 
        subject_filename = os.path.basename(args.subject).split(".")[0]
        database = "database/" + str(subject_filename)
        multifasta2database(args.subject, args.sequence_type, subject_filename)

    blast_compute(args.query, database, args.sequence_type, E_VALUE, fasta=args.fasta, headers=args.headers)
    

if __name__ == '__main__':
    main()