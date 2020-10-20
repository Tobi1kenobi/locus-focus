library(here)
library(tidyverse)
library(EnsDb.Hsapiens.v79)

# Read in results, and get Gene IDs from Ensembl names

pops.16 <- read_tsv(here('results/IBD_DeLange2017.16.results'))
pops.16 <- pops.16 %>% 
  mutate(GENEID = ENSGID) %>% 
  left_join(tibble(ensembldb::select(EnsDb.Hsapiens.v79, keys= pops.16$ENSGID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))), by="GENEID")

pops.20 <- read_tsv(here('results/IBD_DeLange2017.20.results')) 
pop.20 <- pop.20 %>% 
  mutate(GENEID = ENSGID) %>% 
  left_join(tibble(ensembldb::select(EnsDb.Hsapiens.v79, keys= pops.20$ENSGID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))), by="GENEID")

pops.21 <- read_tsv(here('results/IBD_DeLange2017.21.results'))
pops.21 <-  pops.21 %>% 
  mutate(GENEID = ENSGID) %>% 
  left_join(tibble(ensembldb::select(EnsDb.Hsapiens.v79, keys= pops.21$ENSGID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))), by="GENEID")

pops.22 <- read_tsv(here('results/IBD_DeLange2017.22.results'))
pops.22 <- pops.22 %>% 
  mutate(GENEID = ENSGID) %>% 
  left_join(tibble(ensembldb::select(EnsDb.Hsapiens.v79, keys= pops.22$ENSGID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))), by="GENEID")
