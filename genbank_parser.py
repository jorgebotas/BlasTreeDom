#!/usr/bin/env python

import os
import sys
import csv

from Bio import SeqIO

import system

def genBank2multifasta(genBank, sequence_type, output_file = None):
    """ Parse GenBank file and extract all CDS.
        Return or create file containing nucleotide or protein sequences in multifasta depending on sequence_type """
    multifasta = ""
    with open(genBank, 'r') as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            seq = record.seq # DNA sequence
            for feature in record.features:
                start = feature.location.start # Start position 
                end = feature.location.end # End position
                strand = feature.location.strand # Positive or negative strand 
                if feature.type == 'CDS':
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
    # Create .fasta file containing nucleotide/protein sequences
    if output_file:
        output = open(str(output_file) + ".fasta", 'r+')
        output.write(multifasta)
        output.close()
        return
    # If no output_file is inputed, return multifasta 
    return multifasta


def gb_dir2multifasta(genBank_dir, sequence_type, output_filename):
    """ Generate multifasta with all sequences in all genBank files in input directory """
    genBank_list = system.list_files(genBank_dir)
    genBank_multifasta = ""
    for genBank_doc in genBank_list:
        genBank_multifasta += genBank2multifasta(genBank_doc, sequence_type)
    output_file = open(output_filename, "w")
    output_file.write(genBank_multifasta)
    output_file.close()
    return


def field_parser(input_gb, fields = ['locus_tag', 'gene', 'protein_id', 'EC_number', 'product', 'db_xref'], output_file=None):
    """ Generate tsv file containing input fields from CDS in input genbank file """
    if not output_file:
        output_file = str(os.path.basename(input_gb)) + "_fields.tsv"
    tsv_file = open(output_file, 'w')
    tsv_writer = csv.writer(tsv_file, delimiter = '\t')
    if fields:
        tsv_writer.writerow(fields)
    with open(input_gb, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    data = []
                    if not fields:
                        fields = feature.qualifiers.keys()
                    for field in fields:
                        try:
                            value = str(feature.qualifiers[field])
                            for ch in "[]'":
                                value = value.strip(ch)
                            data.append(value)
                        except:
                            data.append("NA")
                    tsv_writer.writerow(data)
    return


def main():
    input_gb = sys.argv[1]
    gb_basename = os.path.basename(input_gb)
    # sequence_type = sys.argv[2]
    # genBank2multifasta(input_gb, sequence_type)

    # Default fields to parse
    fields = ['locus_tag', 'gene', 'protein_id', 'EC_number', 'product', 'db_xref']
    field_parser(input_gb, fields, str(gb_basename) + "_fields.tsv")

    sys.exit()


if __name__ == '__main__':
    main()