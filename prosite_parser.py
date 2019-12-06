#!/usr/bin/env python

import re

from Bio.ExPASy import Prosite,Prodoc


def dat_parser(sequence, fields = ["name", "accession", "description", "pattern"]):
    """ Finds domain hits from prosite.dat in input sequence """
    hits = []
    handle = open("prosite_docs/prosite.dat", "r")
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
    # con este script podemos parsear el archivo .doc
    handle = open("prosite.doc")
    records = Prodoc.parse(handle)
    for record in records:
        print('accession' + str(record.accession))
        print('refs ' + str(record.prosite_refs))
        # print(record.text)
        #print(record.references)
    return

def main():
    fasta = open("data/PBP1_staphylococcus.fasta", "r")
    sequence = fasta.read()[64:].strip('\n')
    print(dat_parser(sequence))



if __name__ == '__main__':
    main()