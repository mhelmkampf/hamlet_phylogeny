### ============================================================================
### Hamlet phylogeny: R code to plot Fig. S2 (ASTRAL summary tree)
### Created Oct 2024 by Martin Helmkampf
### ============================================================================


### Preparations
library(ape)
library(treeio)
library(phangorn)
library(ggtree)
library(tidyverse)

# Set path to root directory of git repository "hamlet_phylogeny"

dev.off()
rm(list=ls())


### Read in tree
phylo <- read.tree("results/trees/phylo2k_s27_astral.tre")


### Set terminal branches to constant length
phylo$edge.length <- replace(phylo$edge.length, phylo$edge.length == "NaN", 0.025)


### Root tree
outgroup <- c("tab", "tig", "tor")
(outsamples <- grep(paste(outgroup, collapse = "|"), phylo$tip.label, value = TRUE))   # identify samples matching outgroup pattern

rooted <- root(phy = phylo, outgroup = outsamples, edgelabel = TRUE)


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
longg <- which(tree_data$branch.length > median(tree_data$branch.length) * 5)


### Compress long branches (requires rooted phylo object index)
longr <- which(rooted$edge.length > median(rooted$edge.length) * 5)
rooted$edge.length[longr] <- rooted$edge.length[longr] * 0.05


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
           support = as.numeric(if_else(!isTip, label, "NA")) * 100,
           support_class = cut(as.numeric(support), c(0, 50, 70, 90, 100)) %>%
             as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")),
           branch_type = case_when(node %in% longg ~ "broken", 
                                   TRUE ~ "whole")
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


### Draw tree
(t <- ggtree(tr = tree, 
             layout = 'rectangular',
             size = 0.4,
             color = scol,
             aes(linetype = branch_type)) +
    geom_tiplab(color = "gray60", 
                size = 1.5, 
                hjust = -0.1,
                offset = 0.0005) +
    geom_tippoint(aes(color = region), 
                  size = 1,
                  alpha = 0.5) +
    geom_nodepoint(data = tree %>% filter(!isTip, support_class != "(0,50]"),
                   aes(fill = support_class,
                       size = support_class),
                   shape = 21,
                   color = "gray20") +
    # labs(title = "Hamlet phylogeny, genomic windows summary",
    #      subtitle = "ASTRAL, 2000 random 5 kb windows, local trees inferred with IQ-TREE, default settings") +
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
    ggtree::geom_treescale(width = 0.1,
                           offset = -4, 
                           fontsize = 3, 
                           color = "gray60",
                           x = 0, y = 50) +
    theme(legend.position = c(0.1, 0.85),
          legend.text = element_text(color = "gray20"),
          legend.title = element_text(color = "gray20"),
          legend.spacing.y = unit(5, "mm")
    )
)


### Save to file
ggsave(plot = t,
       filename = "figures/FigS2_astral.pdf",
       width = 9,
       height = 14,
       device = cairo_pdf,
       bg = "transparent")



### ============================================================================
### Alternative figures versions

# (t2 <- ggtree(tr = tree, 
#              layout = 'rectangular',
#              size = 0.5,
#              color = scol,
#              aes(linetype = branch_type)) +
#    geom_nodepoint(data = tree %>% filter(!isTip, support_class != "(0,50]"),
#                   aes(fill = support_class,
#                       size = support_class),
#                   shape = 21,
#                   color = "gray20") +
#    # geom_text(data = tree %>% filter(!isTip, support >= 50), aes(label = support), color = "gray40", size = 2, hjust = -0.25) +
#    scale_color_manual(values = c("olivedrab4", "royalblue1", "coral2", "gray60")) +
#    scale_fill_manual(values = c(`(0,50]`   = "transparent",
#                                 `(50,70]`  = "white",
#                                 `(70,90]`  = "gray",
#                                 `(90,100]` = "black"),
#                      drop = FALSE) +
#    scale_size_manual(values = c(`(0,50]`   = 0,
#                                 `(50,70]`  = 1,
#                                 `(70,90]`  = 1,
#                                 `(90,100]` = 1),
#                      na.value = 0,
#                      drop = FALSE) +
#    scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
#    guides(color = guide_legend(title = "Region"),
#           fill = guide_legend(title = "Support", title.position = "top"),
#           size = guide_legend(title = "Support", title.position = "top")
#    ) +
#    ggtree::geom_treescale(width = 0.1,
#                           offset = -4, 
#                           fontsize = 3, 
#                           color = "gray60",
#                           x = 0, y = 50) +
#    theme(legend.position = c(0.12, 0.3),
#          plot.title = element_text(size = 14, color = "gray20", face = "bold", vjust = 3, margin = margin(t = 10)),
#          plot.subtitle = element_text(size = 10, color = "gray20", margin = margin(t = 5, b = 15)),
#          legend.text = element_text(color = "gray20"),
#          legend.title = element_text(color = "gray20"),
#          legend.spacing.y = unit(5, "mm")
#    )
# )
# 
# 
# ### Save to file
# ggsave(plot = t2,
#        filename = "FigS2_astral2.pdf",
#        width = 9,
#        height = 14,
#        device = cairo_pdf,
#        bg = "transparent")
