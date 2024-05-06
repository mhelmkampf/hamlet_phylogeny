### ============================================================================
### Hamlet phylogeny: R code to plot Fig. S6 (Admixture ancestry proportions and
### GNN proportions)
### Created Nov 2023 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)
library(cowplot)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Read in admixture proportions
admix2 <- read_delim("results/admix/AdmcvProp_phyps2e_m2n1l5_k2.tsv", delim = " ", col_names = "Sample") %>%
    mutate(Species = str_sub(Sample, -6, -4), 
           Location = str_sub(Sample, -3, -1)
           )


### Set species colors
colstandard <- read_tsv("metadata/species_colors.tsv",
                        col_types = "cc",
                        col_names = TRUE)


### Define regions
gulf <- c("liz", "tam", "ala", "arc", "are", "flk")
west_carib <- c("bel", "boc", "gun", "hon", "san", "qui")
east_carib <- c("bar", "hai", "pri")


### Pivot to long format
long2 <- pivot_longer(admix2, cols = starts_with("X"), 
                      names_to = "Ancestry", 
                      values_to = "Proportion")


### Order samples by ancestry proportions
order <- long2 %>%
  filter(Ancestry == "X3") %>%
  arrange(desc(Proportion)) %>%
  pull(Sample)


### Plot admixture proportions
(pa <- long2 %>%
    mutate(Sample = factor(Sample, levels = unique(Sample))) %>%
    ggplot(aes(x = fct_relevel(Sample, order), y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity", alpha = 0.9) +
    labs(x = NULL,
         y = "Ancestry proportions") +
    scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
          axis.text.x = element_blank(),
          # axis.text.x = element_text(hjust = 1, vjust = 1, angle = 35, size = 1,
          #                            margin = unit(c(-2, 0, 0, 0), "mm")),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)


### Read in GNN proportions
(gnn <- read_delim("results/gnn/gnn_phylo2e_LG02_cld_n1.csv", delim = ",", 
                   col_names = TRUE) %>%
    rename(Clade = clade,
           GNN_small = small,
           GNN_large = large) %>%
    mutate(Haplotype = paste0(Individual, "-", `Sample node`),
           Species = str_sub(Individual, -6, -4),
           Location = str_sub(Individual, -3, -1),
           Population = str_sub(Individual, -6, -1)
    ) %>%
    arrange(desc(GNN_small)) %>%
    select(Haplotype, Population, Clade, GNN_small, GNN_large)
)


### Average across haplotypes per sample
per_sample <- gnn %>%
  mutate(Sample = substr(Haplotype, 1, nchar(Haplotype) - 4)) %>%
  group_by(Sample) %>%
  mutate(Avg_GNN_small = mean(GNN_small),
         Avg_GNN_large = mean(GNN_large)) %>%
  ungroup() %>%
  distinct(Sample, .keep_all = TRUE)  # keep only one row per Haplotype


### Pivot to long format
longg <- pivot_longer(per_sample,
                      cols = starts_with("Avg"),
                      names_to = "GNN_clade",
                      values_to = "Proportion") %>%
  select(-c("Haplotype", "GNN_small", "GNN_large"))


### Plot GNN proportions (ordered as by Admixture ancestry proportions)
(pg <- longg %>%
    arrange(GNN_clade, Proportion) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample))) %>%
    ggplot(aes(x = fct_relevel(Sample, order), y = Proportion, fill = GNN_clade)) +
    geom_bar(position = "fill", stat = "identity", alpha = 0.9) +
    labs(x = NULL,
         y = "GNN proportions") +
    scale_fill_manual(values = c("peru", "#837ABE")) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
          axis.text.x = element_blank(),
          # axis.text.x = element_text(hjust = 1, vjust = 1, angle = 35, size = 1,
          #                            margin = unit(c(-2, 0, 0, 0), "mm")),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)


### Create species color table
scol <- as_tibble(order) %>%
  rename(Sample = value) %>%
  mutate(Position = seq(1:length(order)),
         Species = str_sub(Sample, -6, -4)) %>%
  left_join(., colstandard, by = join_by(Species))


### Plot species colors
(s <- ggplot(scol, aes(xmin = Position, xmax = Position + 1, 
                       ymin = -Inf, ymax = Inf)) +
    geom_rect(aes(fill = Species)) +
    scale_x_discrete(breaks = seq(1:length(order))) +
    scale_y_continuous(breaks = .5, limits = c(0, 1)) +
    scale_fill_manual(values = colstandard$Color) +
    labs(x = NULL, 
         y = "Species") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5,
                                      color = "gray20", size = 9,
                                      margin = margin(r = 0)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(3, 0, 3, 0), "mm")
    )
)


### Create species color table
rcol <- as_tibble(order) %>%
  rename(Sample = value) %>%
  mutate(Position = seq(1:length(order)),
         Location = str_sub(Sample, -3, -1),
         Region = case_when(
           Location %in% gulf ~ "Gulf of Mexico",
           Location %in% west_carib ~ "Western Caribbean",
           Location %in% east_carib ~ "Eastern Caribbean")
         )


### Plot region colors
(r <- ggplot(rcol, aes(xmin = Position, xmax = Position + 1, 
                       ymin = -Inf, ymax = Inf)) +
    geom_rect(aes(fill = Region)) +
    scale_x_discrete(breaks = seq(1:length(order))) +
    scale_y_continuous(breaks = .5, limits = c(0, 1)) +
    scale_fill_manual(breaks = c("Gulf of Mexico", "Western Caribbean", "Eastern Caribbean"),
                      values = c("coral2", "royalblue1", "olivedrab4")) +
    labs(x = NULL, 
         y = "Region") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5,
                                      color = "gray20", size = 9,
                                      margin = margin(r = 0)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(3, 0, 3, 0), "mm")
    )
)


### Combine panels
(p_comb <- plot_grid(pa, pg, s, r,
                     align = "v", 
                     nrow = 4, 
                     rel_heights = c(1, 1, 0.3, 0.3),
                     scale = c(1, 1, 1.002, 1.002)) +
    theme(plot.margin = margin(b = 5))
)


### Save as PDF
ggsave(
  filename = "figures/FigS6_admixGnn.pdf",
  plot = p_comb,
  width = 20,
  height = 11,
  units = "cm",
  device = cairo_pdf
)



### ============================================================================
### Alternative figure versions

# library(patchwork)
# 
# (p_comb <- pa / pg +
#   plot_annotation(tag_levels = c("a")) +
#   theme(plot.tag.position = c(0, 1),
#         plot.tag = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
# )
