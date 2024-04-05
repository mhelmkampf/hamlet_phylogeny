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

# Add Betancur metadata to table, find missing Ross & 06 accessions


### Identify samples to submit to ENA
to_submit <- samples %>%
  filter(is.na(samples$Accession))


### ENA Webin account
# Webin-67244


### Check by set
(sets <- read_tsv("~/Documents/Projects/1_Resources/0_Metadata/seqdata_wgs18.tsv", col_names = TRUE))

main <- samples %>%
  left_join(., sets, by = join_by(SampleID), keep = FALSE)

(no_acc <- main %>%
  filter(is.na(main$Accession)) %>%
  select(SampleID, Sample, Dataset) %>%
  arrange(Dataset)
)

filter(no_acc, Dataset != "phylo2")

  