#!/usr/bin/env python
import os
import re
from subprocess import call, PIPE, Popen

from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd



def list_all(input_list):
    """ Recursively list all the files in given list of files and directories """
    expanded_list = []
    regex = re.compile(r'^\.') # Remove hidden files when listing a directory
    for item in input_list:
        if os.path.isdir(item):
            ls_dir = [item.rstrip('/')+'/'+file for file in os.listdir(item) if not regex.match(file)]
            expanded_list.extend(list_all(ls_dir))
        else:
            expanded_list.append(item)
    return expanded_list


def merge_files(input_files, output_dir, output_filename):
    """ Merge files from input list of files and/or directories """
    with open(output_dir.rstrip('/')+'/'+os.path.basename(output_filename), 'w') as output_file:
        for file in list_all(input_files):
            with open(file, 'r') as text:
                output_file.write(text.read())
    return 

def fasta2dirs(fasta_file, output_dir):
    """ Generate directory inside output_dir for each sequence id in FASTA file """
    if not os.path.isdir(output_dir): os.mkdir(output_dir)
    with open(fasta_file, 'r') as fasta:
        for title, dummy_seq in SimpleFastaParser(fasta):
            identifier = title.split(None, 1)[0]
            directory = output_dir.rstrip('/')+'/'+identifier
            if not os.path.isdir(directory): os.mkdir(directory)
    return

def fasta2tsv(fasta_file, output_dir, output_filename, fields=['id', 'seq'], seq_length=False):
    """ Create dataframe from FASTA file """
    ids = []
    seqs = []
    with open(fasta_file, 'r') as fasta:
        for title, sequence in SimpleFastaParser(fasta):
            ids.append(title.split(None, 1)[0])
            seqs.append(sequence)
    identifiers = pd.Series(ids, name='id')
    sequences = pd.Series(seqs, name='seq')
    df = pd.DataFrame(dict(id=identifiers, seq=sequences))
    fields0=['id', 'seq']
    if seq_length:
        df['length'] = df.seq.apply(len)
        fields0.append('length')
    df.rename(columns=dict(zip(fields0, fields)), inplace=True)
    if not os.path.isdir(output_dir): os.mkdir(output_dir) # Create output directory
    df.to_csv(output_dir.rstrip('/')+'/'+os.path.basename(output_filename), index=False, sep='\t')
    return


def tsv2fasta(tsv_file, output_dir, separate_dirs=False):
    """ Create FASTA file with subject sequences for each query id in tsv_file (dataframe).
        FASTA filename(s) in output_dir or query-specific directory:  unaligned.fasta """
    df = pd.read_csv(tsv_file, delimiter='\t')
    for qseqid in pd.unique(df.qseqid):
        if separate_dirs: 
            filename = output_dir.rstrip('/')+'/{}/unaligned.fasta'.format(qseqid, qseqid)
            os.makedirs(output_dir.rstrip('/')+'/'+qseqid, exist_ok=True)
        else: filename = output_dir.rstrip('/')+'/'+'unaligned.fasta' 
        data = df[df.qseqid == qseqid].reset_index(drop=True)
        with open(filename, 'w') as fasta:
            for idx in range(len(data.sseqid)):
                fasta.write('>'+data.sseqid[idx]+'\n')
                fasta.write(data.sseq[idx]+'\n')
    return


def merge_results(blast_output, genBank_tsv, output_dir, output_filename, include_gb=False):
    """ Merge blast_output, genBank parsed fields and extracted ProSite domain names to create an integrated ouput file """
    blast = pd.read_csv(output_dir.rstrip('/')+'/'+os.path.basename(blast_output), delimiter='\t').drop(columns=['qseqlen'])
    domains = []
    for qid in pd.unique(blast.qseqid):
        data = blast[blast.qseqid == qid]
        domain_file = pd.read_csv(output_dir.rstrip('/')+'/'+qid+'/'+'domains/_domains.tsv', delimiter='\t')
        for sseqid in pd.unique(data.sseqid):
            domains.append(str(pd.unique(domain_file[domain_file.id == sseqid].name)).lstrip('[').rstrip(']'))
    blast.insert(loc=7, column='domains',value=pd.Series(domains, name='domains'), allow_duplicates=True)
    if include_gb:
        fields0 = ['protein_id', 'gene', 'locus_tag', 'EC_number', 'product', 'db_xref']
        fields = ['sseqid', 'gene', 'locus_tag', 'EC_number', 'product', 'db_xref']
        genBank = pd.read_csv(output_dir.rstrip('/')+'/'+os.path.basename(genBank_tsv), delimiter='\t').rename(columns=dict(zip(fields0, fields)))
        merged_df = pd.merge(left=genBank, right=blast, on='sseqid')
        columns = merged_df.columns.values
        columns = np.delete(columns, np.argwhere(columns=='qseqid'))
        columns = np.append(np.array(['qseqid']), columns)
        merged_df = merged_df[columns]
        merged_df.to_csv(output_dir.rstrip('/')+'/'+os.path.basename(output_filename), index=False, sep='\t')
    else:
        merged_df = blast
    merged_df = merged_df.sort_values(by=['qseqid', 'qcovs', 'pident'], ascending=[1, 0, 0])
    merged_df.to_csv(output_dir.rstrip('/')+'/'+os.path.basename(output_filename), index=False, sep='\t')
    return
