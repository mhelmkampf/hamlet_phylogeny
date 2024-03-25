### ============================================================================
### Hamlet phylogeny: R code to plot GNN proportions (tskit)
### Created Mar 2024 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Read in data
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


### Pivot to long format
long <- pivot_longer(gnn,
                     cols = starts_with("GNN"),
                     names_to = "GNN_clade",
                     values_to = "Proportion")


### Plot GNN proportions
(g <- long %>%
  arrange(GNN_clade, Proportion) %>%
  mutate(Haplotype = factor(Haplotype, levels = unique(Haplotype))) %>%
  ggplot(aes(x = Haplotype, y = Proportion, fill = GNN_clade)) +
    geom_bar(position = "fill", stat = "identity", alpha = 0.9) +
    labs(x = "Haplotypes",
         y = "GNN proportions") +
    scale_fill_manual(values = c("peru", "#837ABE")) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
          axis.text.x = element_text(hjust = 1, vjust = 1, angle = 35, size = 1,
                                     margin = unit(c(-2, 0, 0, 0), "mm")),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)


### Save as PDF
ggsave(
  filename = "figures/FigS6b_gnn.pdf",
  plot = g,
  width = 25,
  height = 5,
  units = "cm",
  device = cairo_pdf
)
