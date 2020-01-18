#!/usr/bin/env python

import os
import re

from Bio.ExPASy import Prosite,Prodoc
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
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



def store_domain_info(input_sequence, output_filename, fields=['name', 'accession', 'description', 'pattern'], location=False):
    """ Given an input protein sequence, find ProSite domains and store in output_filename.
        Return domains found """
    output_file = open(output_filename, 'w')
    domains = dat_parser(input_sequence, fields=fields)
    located_domains = []
    output_file.write(str(len(domains))+' domains found.\n\n\n')
    for domain in domains:
        if location:
            matches = re.finditer(domain[3], input_sequence)
            for match in matches:
                located_domains.append(domain + [match.start(), match.end(), match.start() + (match.end() - match.start())/2])
        output_file.writelines([str(field)+'\n' for field in domain]+['\n'])
        text = doc_parser(domain[1])
        output_file.write(text + '\n')
    output_file.close()
    if location: return located_domains
    else: return domains
    ## Create file containing start and end positions of domain in each of the sequences


def extract_domains(input_fasta, output_dir, summary=True):
    """ Given a FASTA file, extract domains of each sequence, store in different files under same directory.
        Create summary file if requested """
    seqids = []
    # Store domain info: 'name', 'accession', 'description', 'pattern', 'start', 'end', 'midpoint'
    columns = ['name', 'accession', 'description', 'pattern', 'start', 'end', 'midpoint']
    domains = [ [] for col in columns ]
    with open(input_fasta, 'r') as fasta:
        for title, sequence in SimpleFastaParser(fasta):
            seqid = title.split(None, 1)[0]
            seq_domains = store_domain_info(input_sequence=sequence, output_filename=output_dir.rstrip('/')+'/'+seqid+'_dominfo.txt', 
                                            fields=columns[:4], location=True)
            for domain in seq_domains:
                for idx in range(len(domains)):
                    domains[idx].append(domain[idx])
            seqids.extend([seqid for dummy in range(len(seq_domains))])
    df = pd.DataFrame()
    df['id'] = pd.Series(seqids, name='id')
    for idx in range(len(domains)):
        df[columns[idx]] = pd.Series(domains[idx], name=columns[idx])
    df.to_csv(output_dir.rstrip('/')+'/'+'_domains.tsv', index=False, sep='\t')



def find_domains(blast_output, output_dir, summary=True):
    """ For each query in blast_output.tsv, extract domains of every sequence in unaligned.fasta """
    df = pd.read_csv(blast_output, delimiter='\t')
    for qid in pd.unique(df.qseqid):
        query_dir = output_dir.rstrip('/')+'/'+qid+'/'
        os.makedirs(query_dir+'domains/', exist_ok=True)
        sseqs_file = query_dir+'/'+'unaligned.fasta'
        ############### EXTRACT DIRECTLY FROM DF ###############
        extract_domains(sseqs_file, query_dir+'domains/', summary=summary)
    return