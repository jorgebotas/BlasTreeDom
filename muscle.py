#!/usr/bin/env python

import os
import sys

from subprocess import call, PIPE, Popen

def multiple_alignment(multifasta, output_filename, query_fasta = None):
    # Perform multiple alignment using MUSCLE
    # If query_fasta provided, include in multiple alignment input
    if query_fasta:
        multifasta_file = open(multifasta, 'r+')
        query = open(query_fasta, 'r')
        multifasta_file.write(query.read())
        multifasta_file.close()
        query.close()
    call(['muscle', '-in', multifasta, '-out', output_filename])
    return


def compute_NJtree(alignment, output_filename = "NJ.tree"):
    # Compute NJ tree using MUSCLE
    # muscle -maketree -in $alignment -out "$directory/NJ_$type.tree" -cluster neighborjoining 1>>$log 2>>$log
    call(['muscle', '-maketree', '-in', alignment, '-out', output_filename, '-cluster', 'neighborjoining'])
    return 
