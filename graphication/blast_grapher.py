#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def blast_plot(blast_output, output_dir, show=False):
    """ Plot blast output for each of the queries provided in input file """
    df = pd.read_csv(blast_output, delimiter='\t')
    df = df.sort_values(['qseqid', 'qcovs', 'pident', 'evalue'], ascending=[1, 0, 0, 0]).reset_index(drop=True)
    df['pident greater than'] = round(df.pident / 10) * 10 # Use pindent to set color accordingly
    queries = pd.unique(df.qseqid)
    for query in queries:
        data = df[df.qseqid == query]
        # with plt.style.context('dark_background'):
        sns.set(context='paper', font_scale=.8, font='times')
        # plt.rc('ytick', labelsize=2)
        fig, ax = plt.subplots(figsize=(5,5), clear=True)
        # Plot totality of query (100% coverage)
        sns.barplot(data=data, x=np.full(len(data.sseqid), 744), y='sseqid' , color="lightgrey") 
        # Plot coverage 
        sns.barplot(data=data, x="qend", y="sseqid", hue='pident greater than', palette=sns.light_palette("green"), dodge=False) #dodge avoids hue shrinkage of width
        # Plot N-terminal if uncovered
        sns.barplot(data=data, x="qstart", y="sseqid", color='lightgrey')
        ax.set(xlabel='Blast overlap', ylabel='subject Accession Number', title='Blast output plot')
        fig.tight_layout()

        fig.savefig("{}{}_blast.png".format(output_dir, query))
        if show:
            plt.show()


# blast_plot("blast_output.tsv")