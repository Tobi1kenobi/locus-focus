import pandas as pd

phrase_list = ['IBD', 'inflammation', 'inflammatory bowel disease', 'Crohns' , 'ulcerative colitis', 'immune', 'gut', 'bowel']
PHRASES = ','.join(phrase_list)
gene_df = pd.read_csv('data/gene_coordinates.GRCh37.ensembl_v91.txt', sep='\t',names=['ENS', 'CHR', 'START', 'END', 'STRAND', 'NAME'])
GENES = list(gene_df['NAME'][1:5])

rule all:
    input:
        "out/gene_distributions.csv"


rule find_counts:
    output:
        "temp/{gene}_phrases.csv"
    conda:
        "conda-env_data-science.yml"
    params:
        phrases = PHRASES
    shell:
        "python pubmed_search.py --gene {wildcards.gene} --phrases '{params.phrases}'"


rule merge_counts:
    input:
        expand("temp/{gene}_phrases.csv",gene=GENES)
    output:
        "out/gene_distributions.csv"
    run:
        shell("echo gene,{phrases} > {{output}}".format(phrases=PHRASES))
        shell("for file in {input}; do cat $file >> {output}; done")
