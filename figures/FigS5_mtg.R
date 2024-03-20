### ============================================================================
### Hamlet phylogeny: R code to plot Fig. S5 (mitochondrial genome phylogeny)
### Created Nov 2023 by Martin Helmkampf
### ============================================================================


### Preparations
library(ape)
library(treeio)
library(ggtree)
library(tidyverse)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list=ls())
dev.off()


### Read in tree in Newick format
phylo <- read.tree("results/trees/phyps2_mtg_GTRG.raxml.support")


### Root tree at midpoint
rooted <- midpoint(phylo)


### Relabel samples
rooted$tip.label <- rooted$tip.label %>%
  str_replace(pattern = "54761nigliz", "54761atlliz") %>%
  str_replace(pattern = "52989nigtam", "52989atltam") %>%
  str_replace(pattern = "27936flotam", "27936atltam") %>%
  str_replace(pattern = "PL17_35puepri", "PL17_35indpri") %>%
  str_replace(pattern = "19294tangun", "19294affgun") %>%
  str_replace(pattern = "62585uniarc", "62585esparc") %>%
  str_replace(pattern = "62571uniarc", "62571esparc") %>%
  str_replace(pattern = "62555uniarc", "62555esparc") %>%
  str_replace(pattern = "27698gutqui", "27698abequi")


### Identify long branches (requires ggtree object index)
tree_data <- ggtree(rooted)$data
longg <- which(tree_data$branch.length > median(tree_data$branch.length) * 20)


### Compress long branches (requires rooted phylo object index)
longr <- which(rooted$edge.length > median(rooted$edge.length) * 20)
rooted$edge.length[longr] <- rooted$edge.length[longr] * 0.1


### Define regions
gulf <- c("liz", "tam", "ala", "arc", "are")
caribbean <- c("bel", "boc", "gun", "hon", "san", "qui")
atlantic <- c("bar", "hai", "pri", "flk")


### Add species / location labels and support categories
(tree <- ggtree(rooted) %>%
    .$data %>%
    mutate(spec = if_else(isTip, str_sub(label, -6, -4), "ungrouped"),
           loc = if_else(isTip, str_sub(label, -3, -1), "ungrouped"),
           region = case_when(
             loc %in% gulf ~ "Gulf of Mexico",
             loc %in% caribbean ~ "Caribbean",
             loc %in% atlantic ~ "Atlantic",
             TRUE ~ "NA"),
           support = as.numeric(if_else(!isTip, label, "NA")),
           support_class = cut(as.numeric(support), c(0, 50, 70, 90, 100)) %>%
             as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")),
           branch_type = case_when(node %in% longg ~ "broken", 
                                   TRUE ~ "whole"),
           mito_nuclear = case_when(label %in% incon ~ "red",
                                    TRUE ~ "gray60")
    )
)


### Set species colors
colstandard <- read_tsv("metadata/species_colors.tsv",
                        col_types = "cc", 
                        col_names = TRUE)

scol <- tree %>%
  rename(Species = spec) %>%
  left_join(colstandard) %>%
  replace_na(list(Color = "gray20")) %>%   # internal branch color
  pull(Color)


### Plot tree
(t <- ggtree(tree, 
             layout = "rectangular",
             size = 0.5,
             color = scol,
             aes(linetype = branch_type)) +
    geom_tiplab(color = "gray60",
                size = 1.5,
                hjust = -0.1) +
    geom_tippoint(aes(color = region), 
                  size = 1,
                  alpha = 1) +
    geom_nodepoint(data = tree %>% filter(!isTip, support_class != "(0,50]"),
                   aes(fill = support_class,
                       size = support_class),
                   shape = 21,
                   color = "gray20") +
    # labs(title = "Hamlet phylogeny, whole mitochondrial genome",
    #      subtitle = "RAxML-NG, thorough search, GTR+G model, 200 nonparametric bootstrap replicates, midpoint root") +
    scale_color_manual(values = c("olivedrab4", "royalblue1", "coral2", "gray60")) +
    scale_fill_manual(values = c(`(0,50]`   = "transparent",
                                 `(50,70]`  = "white",
                                 `(70,90]`  = "gray",
                                 `(90,100]` = "black"),
                      drop = FALSE) +
    scale_size_manual(values = c(`(0,50]`   = 0,
                                 `(50,70]`  = 1,
                                 `(70,90]`  = 1,
                                 `(90,100]` = 1),
                      na.value = 0,
                      drop = FALSE) +
    scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
    guides(color = guide_legend(title = "Region"),
           fill = guide_legend(title = "Support", title.position = "top"),
           size = guide_legend(title = "Support", title.position = "top")
    ) +
    ggtree::geom_treescale(width = 0.001,
                           offset = -4, 
                           fontsize = 3, 
                           color = "gray60",
                           x = 0, y = 255) +
    theme(plot.margin = margin(t = 5, b = 5, unit = "pt"),
          legend.position = c(0.1, 0.91),
          legend.text = element_text(color = "gray20"),
          legend.title = element_text(color = "gray20"),
          legend.spacing.y = unit(5, "mm")
    )
)


### Highlight mito-nuclear incongruencies 
incon <- tree %>%
  filter(label %in% c("PL17_141uniflk", "54650casliz", "62560floarc", "62559floarc")) %>%
  select(label, x, y)


(a <- t +
    geom_segment(data = incon,
                 color = "red",
                 lwd = 1,
                 arrow = arrow(length = unit(0.15, "cm")),
                 aes(x = x + 0.00075,
                     y = y,
                     xend = x + 0.0005,
                     yend = y),
                 inherit.aes = FALSE
                 )
)


### Save plot to file
ggsave(plot = a,
       filename = "figures/FigS5_mtg.pdf",
       width = 10,
       height = 14,
       device = cairo_pdf,
       bg = "transparent")
