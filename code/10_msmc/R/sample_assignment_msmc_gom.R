#!/usr/bin/env Rscript

### ===============================================================================
### R script to create sample groups for MSMC analysis following Hench et al. 2022
### edited by Martin Helmkampf to accommodate phylo2 dataset (pops with n > 2 only)
### relies on functions By Kosmas Hench
### ===============================================================================

library(tidyverse)

# To run from bash command line, execute
# Rscript --vanilla R/sample_assignment_msmc_mh.R R/distribute_samples_msmc_and_cc.R R/cross_cc.R R/sample_info_phylo2e.txt msmc

# Expected working directory: msmc run folder, with R subdirectory

# Uncomment to run from bash command line:
args = commandArgs(trailingOnly=FALSE)
args = args[7:8]

# Uncomment to run from RStudio console:
args <- c("R/distribute_samples_msmc_and_cc.R", "R/cross_cc.R")

distr_script <- as.character(args[1])
cross_script <- as.character(args[2])

print(args)


### Load functions
source(distr_script)
source(cross_script)


### Set seed
set.seed(27678)


### Read in sample list
samples <- read_tsv(file = "../../metadata/ids_phylo2e.txt", col_names = "label") %>%
  mutate(spec = str_sub(label, -6, -4),
         geo = str_sub(label, -3, -1)) %>% 
  filter(!(spec %in% c("tor","tab","tig"))
  )

group_layout <- tibble(n = 1:12,
                       grps_msmc_n = ceiling(n / 4),
                       n_4 = c(0,0,0,1,0,0,1,2,1,2,2,3),
                       n_3 = c(0,0,1,0,1,2,1,0,1,0,1,0),
                       n_2 = c(0,1,0,0,1,0,0,0,1,1,0,0),
                       n_1 = c(1,0,0,0,0,0,0,0,0,0,0,0)
                       )


group_sizes <- samples  %>%
  group_by(spec, geo) %>%
  summarise(n = length(geo),
            inds = list(label)) %>%
  left_join(., group_layout) %>% 
  ungroup() %>%
  filter(n >= 3) %>%   # remove populations with n < 3
  select(-n_2, -n_1)   # remove n_2 and n_1 categories


msmc_grouping <- group_sizes %>%
  select(inds, n_4, n_3, spec, geo) %>%
  pmap(collapse_samples_msmcs) %>% 
  bind_rows() %>%
  mutate(msmc_run = row_number()) %>%
  select(msmc_run, spec:group_nr, group_size, samples)


write_delim(x = msmc_grouping, file = str_c("Ne_grouping_phylo2e_n3.tsv"), delim = "\t")



### ===============================================================================
### CC analysis small vs large clade

gulf <- c("liz", "tam", "ala", "arc", "are", "flk")

samples_gom <- samples %>%
  filter(geo %in% gulf)  # filter to Gulf only

samples_gom$geo <- "gom"
  

group_sizes_cc <- samples_gom  %>%
  group_by(spec, geo) %>%
  summarise(n = length(geo),
            inds = list(label)) %>%
  left_join(., group_layout) %>% 
  ungroup() %>%
  filter(n >= 2) %>%  # remove populations with n < 2
  select(-n_1)


cc_grouping <- group_sizes_cc %>%
  select(inds, spec, geo) %>%
  pmap(collapse_samples_cc) %>%
  bind_rows() %>%
  select(spec:group_nr, samples)


cc_samples <- cc_grouping %>%
  mutate(id = str_c(spec, "_", geo,"_", group_nr)) %>%
  select(id, samples)


cc_tibbles <- cc_grouping %>%
  group_by(spec, geo) %>%
  dplyr::count() %>%
  group_by(geo) %>%
  add_count(name = "n_sp") %>%
  filter(n_sp > 2) %>%   # remove locations with only 1 species
  summarise(content = list(tibble(spec = spec, n = n))) %>%
  bind_cols(., .$content %>% map(paste_groups_cc) %>% tibble) %>%
  setNames(., nm = c("geo", "content", "contrasts"))


cc_output <- cc_tibbles %>%
  select(-content) %>%
  pmap(spread_tibble_cc) %>%
  bind_rows() %>%
  mutate(id1 = str_c(spec_1, "_", geo, "_", group_1),
         id2 = str_c(spec_2, "_", geo, "_", group_2)) %>%
  left_join(., cc_samples %>% set_names(., nm = c("id1","samples_1")))%>%
  left_join(., cc_samples %>% set_names(., nm = c("id2","samples_2"))) %>%
  mutate(run_nr = row_number()) %>%
  select(run_nr, geo, spec_1, spec_2, contrast_nr, samples_1, samples_2)


write_delim(x = cc_output, file = "cc_grouping_phylo2e_gom.tsv", delim = "\t")
