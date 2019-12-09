#!/usr/bin/env python

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def calculate_pident_hue(input_df):
    df = input_df
    df['pident_hue'] = df.pident // 10
    # pident_range = [15, 25, 50, 75, 90, 97, 100]
    # pident_hue = 
    # df['pident_hue'] = pident_hue
    # for pident in df.pident:
    #     for pid in pident_range:
    #         if  pident >= pid:
    #             pident_hue.append(pid)
    #     print(df.loc[df.pident >= pident])
    # df.loc[df.pident < pident_range[-1]]['pident_hue'] = 1 
    return df

def data_handler(input_file):

    df = pd.read_csv(input_file, delimiter='\t')
    df = df.sort_values(['qseqid', 'qcovs', 'pident', 'evalue'], ascending=[1, 0, 0, 0]).reset_index(drop=True)
    queries = pd.unique(df.qseqid)
    for query in queries:
        data = df[df.qseqid == query]
        data = calculate_pident_hue(data)
        query_plot = sns.barplot(data=data, x=np.full(len(data.sseqid), 744), y='sseqid' , color="grey")
        coverage_plot = sns.barplot(data=data, x="qend", y="sseqid", hue='pident_hue', palette=sns.light_palette("green"), dodge=False) #dodge avoids hue shrinkage of width
        n_term = sns.barplot(data=data, x="qstart", y="sseqid", color='grey')


    # plt.savefig("/Users/blackhoodie/Desktop/hola.png")
    plt.show()



data_handler("blast_output.tsv")