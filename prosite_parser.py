#!/usr/bin/env python

import re

from Bio.ExPASy import Prosite,Prodoc
import pandas as pd


def dat_parser(sequence, fields=["name", "accession", "description", "pattern"]):
    """ Finds domain hits from prosite.dat in input sequence """
    hits = []
    with open("prosite_files/prosite.dat", "r") as handle:
        records = Prosite.parse(handle)
        for record in records:
            pattern = record.pattern.strip('.')
            pattern_replacements = {'-' : '',
                                    '{' : '[^', # {X} = [^X]
                                    '}' : ']',
                                    '(' : '{', # (from, to) = {from, to}
                                    ')' : '}',
                                    'X' : '.', # x, X = any (.)
                                    'x' : '.',
                                    '<' : '^', # < = N-terminal
                                    '>' : '$' # > = C-terminal
                                    }
            # Transform ProSite patterns to regular expressions readable by re module
            for pat, repl in pattern_replacements.items():
                pattern = pattern.replace(pat, repl)
            if pattern != "" and re.search(pattern, sequence):
                hits.append([record.name, record.pdoc, record.description, pattern]) # MODIFY
    return hits


def doc_parser(accession):
    """ Returns information on domain with input accession number """
    with open("prosite_files/prosite.doc") as handle:
        records = Prodoc.parse(handle)
        for record in records:
            if accession == record.accession:
                return record.text



def store_domain_info(input_sequence, output_filename, fields=["name", "accession", "description", "pattern"]):
    """ Given an input protein sequence, find ProSite domains and store in output_filename.
    Return num_domains found """
    output_file = open(output_filename, 'w')
    domains = dat_parser(input_sequence, fields=fields)
    output_file.write(str(len(domains))+' domains found.\n\n\n')
    for domain in domains:
        output_file.writelines([field+'\n' for field in domain]+['\n'])
        text = doc_parser(domain[1])
        output_file.write(text + '\n')
    output_file.close()
    return len(domains)
    ## Create file containing start and end positions of domain in each of the sequences



def extract_domains(input_fasta, output_dir, summary=True):
    """ Given a FASTA file, extract domains of each sequence, store in different files under same directory.
        Create summary file if requested """
    total_domains = 0
    n_seq = 0
    with open(input_fasta, 'r') as input_file:
        lines = input_file.readlines()
        for idx in range(len(lines)):
            if lines[idx][0] == '>':
                n_seq += 1
                total_domains += store_domain_info(input_sequence=lines[idx+1].strip('\n'), output_filename=output_dir+lines[idx][1:].strip('\n'))
            else:
                pass
        
    # CREATE SUMMARY FILE

def find_domains(blast_output, output_dir, summary=True):
    """ For each query in blast_output.tsv, extract domains of every sequence in query_sseqs.fasta """
    df = pd.read_csv(blast_output, delimiter='\t')
    for qid in pd.unique(df.qseqid):
        query_dir = output_dir.rstrip('/')+'/'+qid+'/'
        sseqs_file = query_dir+qid+'_sseqs.fasta'
        ############### EXTRACT DIRECTLY FROM DF ###############
        extract_domains(sseqs_file, query_dir, summary=summary)
    return


def main():
    """ Finds ProSite domains in protein and stores pertaining information """
    fasta = open("data/PBP1_staphylococcus.fasta", "r")
    sequence = fasta.read()[64:].strip('\n')
    # print(dat_parser(sequence))
    store_domain_info(sequence, "domain_info.txt")
    # print(doc_parser('PS00001'))
    # re.findall(r"G[^EDRKHPFYW].{2}[STAGCN][^P]", sequence)



if __name__ == '__main__':
    main()