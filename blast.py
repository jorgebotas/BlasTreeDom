#!/usr/bin/env python

import os
import sys

from subprocess import call, PIPE, Popen


def multifasta2database(multifasta, sequence_type, output_filename):
    # Create database from multifasta using makeblastdb
    # Create directory for database files if it does not already exist
    database_path = os.path.abspath("database")
    try:
        os.mkdir(database_path)
    except:
        pass
    call(['makeblastdb', '-in', multifasta, '-dbtype', sequence_type, '-out', 'database/' + str(output_filename)])
    return 


def blast_compute(query_fasta, database_path, sequence_type, e_value,  output_filename = "blast_output.fasta"):
    # Perform blastp or blastn analysis for protein or nucleotide sequences respectively
    if sequence_type == "prot":
        blast_type = 'blastp'
    elif sequence_type == "nucl":
        blast_type = 'blastn'

    blast = Popen([blast_type, '-query', query_fasta, '-db', database_path, '-evalue', e_value, '-outfmt', '6 sseqid sseq evalue pident'], stdout=PIPE, stderr=PIPE) 
    blast_output = blast.stdout.read().decode('utf-8')
    blast.stderr.close()
    blast.stdout.close()

    output_file = open(output_filename, 'w')
    # Generate dictionary with blastp output sequences
    # keys = seq id's 
    # values = sequences
    # blast_dictionary = {}    
    for hit in blast_output.splitlines():
        sseqid = hit.split()[0]
        sseq = hit.split()[1]
        # blast_dictionary[sseqid] = sseq
        output_file.write(">" + str(sseqid) + "\n")
        output_file.write(sseq + "\n")
    output_file.close()
    return # blast_dictionary