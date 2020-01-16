#!/usr/bin/env python

import argparse
import os
import sys

import numpy as np
import pandas as pd
from subprocess import call, PIPE, Popen

import file_handler as fh


def multifasta2database(multifasta, sequence_type, output_dir, output_filename='subject', log='/dev/null'):
    """ Create database from multifasta using makeblastdb """
    # Create directory for database files if it does not already exist
    database_path = output_dir.rstrip('/')+'/database/'
    if not os.path.isdir(database_path): os.mkdir(database_path)
    output = database_path+os.path.basename(output_filename)
    call(['makeblastdb', '-in', multifasta, '-dbtype', sequence_type, '-out', output], stdout=open(log), stderr=open(log)) # -parse_seqids could be added
    return 


def save_multifasta(input_file = "blast_output.tsv", output_filename = "blast_output.fasta"):
    blast_output = pd.read_csv(input_file, delimiter='\t')
    output_file = open(output_filename, 'w')
    for dummy_idx, hit in blast_output.iterrows():
        output_file.write(">" + str(hit.sseqid) + "\n")
        output_file.write(str(hit.sseq) + "\n")
    output_file.close
    return


def blast_compute(query_fasta, database_path, sequence_type, e_value, cov_threshold=0,  pident_threshold=0,
                  outfmt='6 qseqid sseqid qcovs qstart qend pident evalue', 
                  output_filename = "blast_output", log='/dev/null'):
    """ Perform blastp or blastn analysis for protein or nucleotide sequences respectively.
        Output filtered by query coverage, identity percentage and e-value thresholds """
    if sequence_type == "prot":
        blast_type = 'blastp'
    elif sequence_type == "nucl":
        blast_type = 'blastn'

    call([blast_type, '-query', query_fasta, '-db', database_path, '-evalue', e_value, '-out', output_filename + '.tsv', '-outfmt', outfmt], stderr=open(log)) # , '-logfile', log
    #, stderr=open(os.path.abspath(log))
    # os.remove(log+'.perf') # log.perf file created when calling blast with -logfile argument
    # Include headers in output .tsv file
    first_line = outfmt[2:].replace(" ", "\t")
    blast = open(output_filename + '.tsv', 'r')
    blast_output = blast.read()
    blast.close()
    output = open(output_filename + '.tsv', 'w')
    output.write(first_line + '\n')
    output.write(blast_output)
    output.close()
    # Filter output by coverage and identity percentage thresholds
    if pident_threshold == None: pident_threshold = 0
    if cov_threshold == None: cov_threshold = 0 
    df = pd.read_csv(output_filename + '.tsv', delimiter='\t')
    filtered_output = df.loc[(df.pident >= pident_threshold) & (df.qcovs >= cov_threshold)]
    filtered_output.to_csv(output_filename + '.tsv', sep='\t', index=False)


def retrieve_seqs(query_fasta, subject_multifasta, blast_output, output_dir, output_filename, remove_files=False):
    """ Generate multifasta file with complete hit subject sequences for each query in blast_output file """ 
    if not os.path.isdir(output_dir): os.mkdir(output_dir) # Create output directory
    subject_tsv = os.path.basename(subject_multifasta).rsplit('.', 1)[0]+'.tsv'
    query_tsv = os.path.basename(query_fasta).rsplit('.', 1)[0]+'.tsv'
    fh.fasta2tsv(fasta_file=query_fasta, output_dir=output_dir, output_filename=query_tsv, fields=['qseqid', 'qseq', 'qseqlen'], seq_length=True)
    fh.fasta2tsv(fasta_file=subject_multifasta, output_dir=output_dir, output_filename=subject_tsv, fields=['sseqid', 'sseq'])
    blast = pd.read_csv(blast_output, delimiter='\t')
    sseqs = pd.read_csv(output_dir.rstrip('/')+'/'+subject_tsv, delimiter='\t')
    qseqs = pd.read_csv(output_dir.rstrip('/')+'/'+query_tsv, delimiter='\t')
    # Merge dataframes by sseqid
    blast_sseqs = pd.merge(left=blast, right=sseqs, on='sseqid')
    merged = pd.merge(left=blast_sseqs, right=qseqs, on='qseqid')
    merged.to_csv(output_dir.rstrip('/')+'/'+os.path.basename(output_filename), index=False, sep='\t')
    if remove_files:
        os.remove(output_dir.rstrip('/')+'/'+query_tsv)
        os.remove(output_dir.rstrip('/')+'/'+subject_tsv)
    return


def main():
    # e-value threshold 
    E_VALUE = "1e-03" 

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-query', type=str, required=True)
    exclusive = arg_parser.add_mutually_exclusive_group(required=True)
    exclusive.add_argument('-subject', type=str)
    exclusive.add_argument('-database', type=str)
    arg_parser.add_argument('-sequence_type', choices=['prot', 'nucl'], required=True)
    arg_parser.add_argument('-pident')
    arg_parser.add_argument('-cov')
    args = arg_parser.parse_args()

    if args.database: 
        database = args.database
    else: 
        subject_filename = os.path.basename(args.subject).split(".")[0]
        database = "database/" + str(subject_filename)
        multifasta2database(args.subject, args.sequence_type, subject_filename)

    blast_compute(query_fasta=args.query, database_path=database, sequence_type=args.sequence_type,
                  e_value=E_VALUE, cov_threshold=args.cov, pident_threshold=args.pident)
    

if __name__ == '__main__':
    main()
