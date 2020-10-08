library(tidyverse)
library(ggplot2)
library(here)

here::dr_here()
# Some misfirings here here from snakemake pipeline, need to run again for ~20 genes, replace na for now
gene_phrases_df <- read_csv(here::here('results/gene-phrases_search.csv'))

gene_only_df <- read_csv(here::here('results/gene-only_search.csv'))

quants <- apply(gene_phrases_df[2:length(gene_phrases_df)], 2, quantile, probs=c(0.9,0.95, .99))
quants

pubmed_score <- function(quant_mat, gene, dataframe){
  gene_phrases_df[gene_phrases_df$gene==gene,]
  gene_score <- 0
  for (percent in rownames(quant_mat))
  {
    greater_thans <- gene_phrases_df[gene_phrases_df$gene==gene,2:length(gene_phrases_df)] > quant_mat[percent,]
    gene_score <- gene_score + sum(greater_thans)
    print(greater_thans)
  }
  return(gene_score)
}

pubmed_score(quant_mat = quants, gene = "NOD2", dataframe = gene_phrases_df)
