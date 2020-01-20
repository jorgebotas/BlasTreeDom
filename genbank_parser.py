#!/usr/bin/env python

import os
import sys
import csv

from Bio import SeqIO
import numpy as np
import pandas as pd

import file_handler as fh

def gb_parser(genBank, sequence_type, output_dir=None, output_file = None, parse_fields=False):
    """ Parse GenBank file and extract all CDS.
        Return or create file containing nucleotide or protein sequences in multifasta depending on sequence_type """
    multifasta = ""
    if parse_fields:
        fields = ['protein_id', 'gene', 'locus_tag', 'EC_number', 'product', 'db_xref']
        data = [ [] for dummy_field in fields ]
    with open(genBank, 'r') as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            seq = record.seq # DNA sequence
            record_name = record.name
            for feature in record.features:
                start = feature.location.start # Start position in genomic sequence
                end = feature.location.end # End position in genomic sequence
                strand = feature.location.strand # Positive or negative strand 
                if feature.type == 'CDS':
                    try:
                        try:
                            multifasta += ">" + feature.qualifiers['protein_id'][0] + " "
                        except:
                            multifasta += ">"
                        try:
                            multifasta += feature.qualifiers['product'][0] + " "
                        except:
                            pass
                        multifasta += record.name + "\n"
                        if sequence_type == "prot":
                            multifasta += feature.qualifiers['translation'][0] + "\n"
                        else:
                            if strand > 0:
                                multifasta += seq[start : end] + "\n"
                            else:
                                # If strand is negative, the coding sequence is the reverse complementary
                                multifasta += seq[start : end].reverse_complement() + "\n"
                        if parse_fields:
                            for idx in range(len(fields)):
                                try:
                                    field = feature.qualifiers[fields[idx]]
                                    for char in "[]'":
                                        field = str(field).strip(char)
                                    data[idx].append(field)
                                except:
                                    data[idx].append('N/A')
                    except: pass
    # Create tsv file containing parsed info from features.qualifiers
    if parse_fields:
        genBank_tsv = output_dir.rstrip('/')+'/'+'_genBank_info.tsv'
        df = pd.DataFrame()
        df['record_name'] = pd.Series(np.full(len(data[0]), record_name), name='record_name')
        for idx in range(len(fields)):
            df[fields[idx]] = pd.Series(data[idx], name=fields[idx])
        # If file already exists, append newly parsed genBank fields
        if os.path.exists(genBank_tsv):
            previous_tsv = pd.read_csv(genBank_tsv, delimiter='\t')
            df = previous_tsv.append(other=df, ignore_index=True) # Keep order of genBank parsing
        df.to_csv(genBank_tsv, index=False, sep='\t')
    # Create .fasta file containing nucleotide/protein sequences
    if output_file:
        output = open(str(output_file) + ".fasta", 'a+')
        output.write(multifasta)
        output.close()
        return
    # If no output_file is inputed, return multifasta 
    return multifasta


def parse_gbs(genBanks, sequence_type, output_dir, output_filename):
    """ Generate multifasta with all sequences in all genBank files in input directory """
    genBank_list = fh.list_all(genBanks)
    genBank_multifasta = ""
    for genBank_doc in genBank_list:
        genBank_multifasta += gb_parser(genBank=genBank_doc, sequence_type=sequence_type, output_dir=output_dir, parse_fields=True)
    output = output_dir.rstrip('/')+'/'+os.path.basename(output_filename)
    with open(output, "w") as  output_file:
        output_file.write(genBank_multifasta)
    return