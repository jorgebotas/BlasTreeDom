#!/usr/bin/env python

import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import prosite_parser as prop


def domain_mapper(input_sequence):
    """ Return dictionary mapping:
    Prosite protein domain match to start and end position in input sequence """
    domains = prop.dat_parser(input_sequence, fields = ["name", "accession", "description", "pattern"])
    # print(len(domains))
    domain_map = {}
    domain_list = []
    for domain in domains:
        matches = re.finditer(domain[-1], input_sequence)
        counter = 0
        for match in matches:
            counter += 1
            # Append end, start and name of domain hit in descending order 
            domain_list.append((match.end(), match.start(), domain[0]))
            domain_map[domain[0]] = match.start() + (match.end() - match.start())/2
        # print(counter, domain)
    return domain_map
    # return sorted(domain_list, reverse=True)
    

def domain_plot(blast_output, show=False):
    """ Plot ProSite protein domains on blast hits from input file """
    df = pd.read_csv(blast_output, delimiter='\t')
    df = df.sort_values(['qseqid', 'qcovs', 'pident', 'evalue'], ascending=[1, 0, 0, 0]).reset_index(drop=True)
    df['lensseq'] = len(df.sseq)
    # domain_mapper('RNGKVLAEDVERYKLVAVIDKKASANSKKPRHVVDKKETAKKLSTVIDMKPEEIEKRLSQKKAFQIEFGRKGTNLTYQDKLKIEKMNLPGISLLPETERFYPNGNFASHLIGRAQKNPDTGELKGALGVEKIFDSYLSGSKGSLRYIHDIWGYIAPNTKKEKQPKRGDDVHLTIDSNIQVFVEEALDGMVERYQPKDLFAVVMDAKTGEILAYSQRPTFNPETGKDFGKKWANDLYQNTYEPGSTFKSYGLAAAIQEGAFDPDKKYKSGHRDIMGSRISDWNRVGWGEIPMSLGFTYSSNTLMMHLQDLVGADKMKSWYERFGFGKSTKGMFDGEAPGQIGWSNELQQKTSSFGQSTTVTPVQMLQAQSAFFNDGNMLKPWFVNSVENPVSKRQFYKGQKQIAGKPITKDTAEKVEKQLDLVVNSKKSHAANYRIDGYEVEGKTGTAQVAAPNGGGYVKGPNPYFVSFMGDAPKKNPKVIVYAGMSLAQKNDQEAYELGVSKAFKPIMENTLKYLNVGKSKDDTSNAEYSKVPDVEGQDKQKAIDNVSAKSLEPVTIGSGTQIKAQSIKAGNKVLPHSKVLLLTDGDLTMPDMSGWTKEDVIAFENLTNIKVNLKGSGFVSHQSISKGQKLTEKDKIDVEFSSENVDSNSTNNSDSNSDDKKKSDSKTDKDKSD')
    df['domains'] = domain_mapper(df.sseq) # Obtain ProSite domains for each hit sequence
    print(df)
    queries = pd.unique(df.qseqid)
    for query in queries:
        data = df[df.qseqid == query]
        color_palette = sns.color_palette('hls', len(data.sseqid))
        print(color_palette)
        sns.set(context='paper', font_scale=.8, font='times')
        fig, ax = plt.subplots(figsize=(10,10), clear=True)
        # Plot totality of query (100% coverage)
        sns.barplot(data=data, x=np.full(len(data.sseqid), 744), y='sseqid' , color="lightgrey") 
        # Plot coverage 
        sns.barplot(data=data, x="qend", y="sseqid", color='grey')
        # PLOT DOMAINS


        # Plot N-terminal
        sns.barplot(data=data, x="qstart", y="sseqid", color='lightgrey')
        ax.set(xlabel='Blast overlap', ylabel='subject Accession Number', title='Blast output plot')
        fig.tight_layout()

        # fig.savefig("/Users/blackhoodie/Desktop/{}_blast.png".format(query))
        if show:
            plt.show()


# domain_plot("blast_output.tsv")
sequence = 'RNGKVLAEDVERYKLVAVIDKKASANSKKPRHVVDKKETAKKLSTVIDMKPEEIEKRLSQKKAFQIEFGRKGTNLTYQDKLKIEKMNLPGISLLPETERFYPNGNFASHLIGRAQKNPDTGELKGALGVEKIFDSYLSGSKGSLRYIHDIWGYIAPNTKKEKQPKRGDDVHLTIDSNIQVFVEEALDGMVERYQPKDLFAVVMDAKTGEILAYSQRPTFNPETGKDFGKKWANDLYQNTYEPGSTFKSYGLAAAIQEGAFDPDKKYKSGHRDIMGSRISDWNRVGWGEIPMSLGFTYSSNTLMMHLQDLVGADKMKSWYERFGFGKSTKGMFDGEAPGQIGWSNELQQKTSSFGQSTTVTPVQMLQAQSAFFNDGNMLKPWFVNSVENPVSKRQFYKGQKQIAGKPITKDTAEKVEKQLDLVVNSKKSHAANYRIDGYEVEGKTGTAQVAAPNGGGYVKGPNPYFVSFMGDAPKKNPKVIVYAGMSLAQKNDQEAYELGVSKAFKPIMENTLKYLNVGKSKDDTSNAEYSKVPDVEGQDKQKAIDNVSAKSLEPVTIGSGTQIKAQSIKAGNKVLPHSKVLLLTDGDLTMPDMSGWTKEDVIAFENLTNIKVNLKGSGFVSHQSISKGQKLTEKDKIDVEFSSENVDSNSTNNSDSNSDDKKKSDSKTDKDKSD'
# print(domain_mapper(sequence))
domains = domain_mapper(sequence)

print(len(sequence))

print(domains.keys())
print(domains.values())

# df = pd.DataFrame(data=domains, columns=['end', 'start', 'domain_name'])

# print(df)

# color_palette = sns.color_palette('hls', len(df.domain_name))
# sns.set(context='paper', font_scale=.8, font='times')
# fig, ax = plt.subplots(figsize=(15,5), clear=True)
# # Plot totality of query (100% coverage)
# sns.barplot(data=df, x="end", y=np.full(len(df.domain_name), 1), color="lightgrey") 
# # Plot coverage 
# # sns.barplot(data=df, x='end', y=np.full(len(df.domain_name), 1), color='grey', orient="h")

# plt.show()

# # https://matplotlib.org/3.1.1/users/event_handling.html#mouse-enter-and-leave