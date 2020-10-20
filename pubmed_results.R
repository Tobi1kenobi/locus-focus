library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(readxl)
library(here)

here::dr_here()

# Read in the "[gene] AND [phrase]" search hits csv
# TODO: Some misfirings here here from snakemake pipeline, need to run again for ~20 genes, replace na for now
gene_phrases_df <- read_csv(here::here('results/gene-phrases_search.csv')) %>% 
  drop_na()

# Read i nthe "[gene]" search hits csv
gene_only_df <- read_csv(here::here('results/gene-only_search.csv'))

# Divide the number of hits for gene+phrase byt the number of hits for just the gene
# TODO: Maybe consider filtering out low hits e.g. "[gene]" returns 1 hit and
# "[gene] AND [phrase]" returns 1 hit = highest possible score
norm_gene_phrases_df <- left_join(gene_phrases_df, gene_only_df,by='gene') %>% 
  group_by(gene) %>% 
  mutate_at(vars(IBD: bowel), funs(./ BACKGROUND)) %>% 
  replace(is.na(.), 0) %>% 
  select(-BACKGROUND)
  
  
# Decide on the quantiles we think are important, here 90th, 95th and 99th, then calculate them for each phrase
quants <- apply(norm_gene_phrases_df[3:length(norm_gene_phrases_df)], 2, quantile, probs=c(0.9,0.95, .99))
quants



pubmed_score <- function(quant_mat, gene, dataframe){
  # Calculcates a score out of 21 for the gene and its likely involvement
  # in IBD based on pubmed searching
  dataframe[dataframe$gene==gene,]
  gene_score <- 0
  for (percent in rownames(quant_mat))
  {
    greater_thans <- dataframe[dataframe$gene==gene,3:length(dataframe)] > quant_mat[percent,]
    gene_score <- gene_score + sum(greater_thans)
    # print(greater_thans) # Shows breakdown of score, comment out when running on lots of genes e.g. many loci as on line 55
  }
  return(gene_score)
}

# Examples
pubmed_score(quant_mat = quants, gene = "NOD2", dataframe = norm_gene_phrases_df)
pubmed_score(quant_mat = quants, gene = "TTN", dataframe = norm_gene_phrases_df)


my_loci <- read_excel(here::here('data/LocusFocusMasterSpreadsheet.xlsx'),skip = 8) %>%
  filter(Assigned=='Tobi')

# For each locus, get a score as described above for all genes
pubmeds <- lapply(strsplit(my_loci$`Genes in locus`, ','),FUN = function(gene_list){
  scores <- rep(-1, length(gene_list))
  names(scores) <- gene_list
  for (gene in gene_list){
    if (gene == '.'){
      return(scores)
    }
    scores[gene] <-  pubmed_score(quant_mat = quants, gene = gene, dataframe = norm_gene_phrases_df)
  }
  return(scores)
}
)


###### PLOTS ########

which_quant <- function(phrase, score){
  # For a given phrase and normalised score, find out which quantile the score falls in 
  phrase_quants <- quants[,phrase]
  bool_quants <- score > phrase_quants
  for (percent in rev(names(bool_quants))){
    if (bool_quants[percent] == TRUE){
      return(percent)
    }
  }
  return(NA)
  
}

# Format the normalised dataframe for plotting
plotting.tib <- norm_gene_phrases_df %>% 
  select(-IBD) %>% 
  pivot_longer(cols=!gene,names_to = "phrase",values_to = "score") %>% 
  rowwise() %>% 
  mutate(Quantile= which_quant(phrase, score))


# Plot as a beeswarm plot
ggplot(data=plotting.tib, aes(x=phrase, y=score, colour=Quantile)) + 
  geom_quasirandom(alpha=0.3) +
  geom_label(data=subset(plotting.tib, gene == 'CD40'),
            aes(x=phrase,y=score,label=gene),position=position_quasirandom())+
  scale_colour_brewer(palette = "YlOrRd",na.value='grey') +
  theme_classic()
