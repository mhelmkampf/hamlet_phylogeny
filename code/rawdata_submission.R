### ============================================================================
### Hamlet phylogeny: R code to handle sample metadata (including for 
### submission of raw data to ENA)
### Created Mar 2024 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Read in data
(samples <- read_tsv("metadata/samples_table.tsv", col_names = TRUE))


### Identify samples to submit to ENA
to_submit <- samples %>%
  filter(is.na(samples$Accession))

to_submit$Sample

### ENA Webin account
Webin-67244
         