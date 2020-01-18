#!/usr/bin/env python    

import os    
import time


def friendly_user_interfase():
    """ Handle Friendly User Interfase """   
    query = None
    unaligned = None
    genBank = None
    multifasta = None
    cov = None
    pident = None
    e_value = None
    graph = False
    
    yes = ['y', 'Y', 'yes', 'Yes']
    no = ['n', 'N', 'no', 'No']
    os.system('clear')
    print('\n\nWelcome to BlasTreeDom!!!')
    time.sleep(.5)
    print('\n\nEntering Friendly User Interfase Mode...\n\n\n')
    time.sleep(1)
    is_unaligned = input('First things first, do you already have a set of unaligned sequences in FASTA format? [y/n]\n\n>  ')
    if is_unaligned in yes:
        unaligned = input('\n\nGood!, type path or drag and drop FASTA file(s) &/or directory with unaligned sequences\n\n>  ').strip(' ').split(' ')
    elif is_unaligned in no:
        query = input('\n\nWell then, type path or drag and drop query file(s) &/or directory!\n\n>  ').strip(' ').split(' ')
        gb = input('\n\nWe would also need some subject sequences, would do you like to extract them from a genBank file(s) (.gbk, .gbff...)?\n\n>  ')
        if gb in yes:
            genBank = input('\n\nThen, type type path or drag and drop those genBank file(s) &/or directory\n\n>  ').strip(' ').split(' ')
        elif gb in no:
            multifasta = input('\n\nWe\'d need some FASTA file(s) &/or directory containing such sequences... Type path or drag and drop them\n\n>  ').strip(' ').split(' ')
        else: 
            print('\n\nIncorrect input, try again!\n\n')
            exit(1)
    else: 
        print('\n\nIncorrect input, try again!\n\n')
        exit(1)
    plot = input('\n\nWould you want to awesomely graph the results???\n\n>  ')
    if plot in yes:
        print('\n\nCorrect answer! :)\n\n')
        graph = True
    elif plot in no:
        print('\n\nWhat a pitty :(\n\n')
        graph = False
    else:
        print('\n\nIncorrect input, try again!\n\n')
        exit(1)
    print('STARTING process...\n\n')
    time.sleep(.5)
    return query, unaligned, genBank, multifasta, cov, pident, e_value, graph
