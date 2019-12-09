#!/usr/bin/env python

import re

from Bio.ExPASy import Prosite,Prodoc


def dat_parser(sequence, fields = ["name", "accession", "description", "pattern"]):
    """ Finds domain hits from prosite.dat in input sequence """
    hits = []
    handle = open("prosite_files/prosite.dat", "r")
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
            hits.append([record.name, record.accession, record.description, pattern])
    handle.close()
    return hits


def doc_parser(accession):
    """ Returns information on domain with input accession number """
    handle = open("prosite_files/prosite.doc")
    records = Prodoc.parse(handle)
    for record in records:
        if accession == record.accession:
            pass
        # print('refs ' + str(record.prosite_refs))
        # print(record.text)
        # print(record.references)
    return


def main():
    """ Finds ProSite domains in protein and stores pertaining information """
    fasta = open("data/PBP1_staphylococcus.fasta", "r")
    sequence = fasta.read()[64:].strip('\n')
    print(dat_parser(sequence))
    # print(doc_parser('PS00001'))
    re.findall(r"G[^EDRKHPFYW].{2}[STAGCN][^P]", sequence)



if __name__ == '__main__':
    main()