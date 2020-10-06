import pandas as pd
import numpy as np
from Bio import Entrez
import plotnine as plt9

Entrez.email = 'oa3@sanger.ac.uk'

all_loci = pd.read_excel('data/LocusFocusMasterSpreadsheet.xlsx', skiprows=8)
gene_list = pd.read_csv('data/gene_coordinates.GRCh37.ensembl_v91.txt', sep='\t',names=['ENS', 'CHR', 'START', 'END', 'STRAND', 'NAME'])
tobi_loci = all_loci[all_loci.Assigned=='Tobi'].copy()
tobi_loci

def search_pubmed(query):
    '''
    Searches pubmed with a provided query
    '''
    handle = Entrez.esearch(db="pubmed", term=query)
    record = Entrez.read(handle)
    handle.close()
    return(record)

def phrase_gene_distributions(phrase_list, gene_list):
    '''
    For a given phrase list and gene list, finds all possible
    '''
    gene_counts_df = pd.DataFrame (gene_list,columns=['Gene'])
    gene_counts_df.set_index('Gene', inplace=True)
    for phrase in phrase_list:
        print(phrase)
        for gene in gene_list:
            try:
                gene_counts_df.at[gene, phrase] = int(search_pubmed('{} AND {}"'.format(phrase, gene))['Count'])
            except:
                print('Could not search for {}'.format(gene))
    
    return gene_counts_df

keywords = ['IBD', 'inflammation', 'inflammatory bowel disease', 'Crohn\'s',
            'Ulcerative colitis', 'immune', 'gut', 'bowel', '']
gene_distributions = phrase_gene_distributions(keywords, gene_list['NAME'])

gene_distributions.to_csv('gene_phrase_distributions.csv')