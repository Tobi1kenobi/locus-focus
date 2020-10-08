import pandas as pd
import numpy as np
from Bio import Entrez
import plotnine as plt9
import argparse

################
## PARAMS
################

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gene", type=str,
                    help="the gene to search for")
parser.add_argument("-p", "--phrases", type=str,
                    help="the comma separated phrases to search for")
args = parser.parse_args()

gene = args.gene
phrases = args.phrases.split(',')


Entrez.email = 'oa3@sanger.ac.uk'

##################
## FUNCTIONS
##################

def search_pubmed(query):
    '''
    Searches pubmed with a provided query
    '''
    handle = Entrez.esearch(db="pubmed", term=query)
    record = Entrez.read(handle)
    handle.close()
    return(record)

def phrases_gene_counts(phrase_list, gene):
    '''
    For a given phrase list and gene, counts how many hits they have when searched together
    '''
    counts_df = pd.DataFrame(columns=phrase_list, index=[gene])

    for phrase in phrase_list:
        if phrase == 'NONE':
            count = int(search_pubmed('{}'.format(gene))['Count'])
            no_phrase_df = pd.DataFrame(columns=['BACKGROUND'], index=[gene])
            no_phrase_df.at[gene, 'BACKGROUND'] = count
            no_phrase_df.to_csv('temp/{}_no-phrase.csv'.format(gene),header=False)
        else:
            count = int(search_pubmed('{} AND {}'.format(phrase, gene))['Count'])
            counts_df.at[gene, phrase] = count
            gene_phrase_counts_df.to_csv('temp/{}_phrases.csv'.format(gene),header=False)
    return none

#############
## MAIN
#############

phrases_gene_counts(phrases, gene)