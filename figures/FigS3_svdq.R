### ============================================================================
### Hamlet phylogeny: R code to plot Fig. S3 (SVDQuartets tree)
### Created Oct 2024 by Martin Helmkampf
### ============================================================================


### Preparations
library(ape)
library(treeio)
library(ggtree)
library(tidyverse)

# Set path to root directory of git repository "hamlet_phylogeny"

dev.off()
rm(list=ls())


### Read in majority consensus tree
phylo <- read.nexus("results/trees/phylo2e_m2k5_svdq_bs200C_maj.nex")
phylo <- phylo$MajRule


### Root tree
outgroup <- c("tab", "tig", "tor")
(outsamples <- grep(paste(outgroup, collapse = "|"), phylo$tip.label, value = TRUE))   # identify samples matching outgroup pattern

rooted <- root(phy = phylo, outgroup = outsamples, edgelabel = TRUE)


### Add species / location labels and support categories
(tree <- ggtree(rooted) %>%
    .$data %>%
    mutate(spec = if_else(isTip, str_sub(label, -6, -4), "ungrouped"),
           loc = if_else(isTip, str_sub(label, -3, -1), "ungrouped"),
           support = as.numeric(if_else(!isTip, label, "NA")),
           support_class = cut(as.numeric(support), c(0, 50, 70, 90, 100)) %>%
             as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"))
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
             size = 0.5,
             color = scol) +
    geom_tiplab(color = "gray60", 
                size = 1.5, 
                hjust = -0.1,
                offset = 0.0005) +
    # geom_text(data = tree %>% filter(!isTip, support >= 50), aes(label = support), color = "gray40", size = 2, hjust = -0.25) +
    geom_nodepoint(data = tree %>% filter(!isTip, support_class != "(0,50]"),
                   aes(fill = support_class,
                       size = support_class),
                   shape = 21,
                   color = "gray20") +
    # labs(title = "Hamlet phylogeny, multispecies coalescent model",
    #      subtitle = "SVDQuartets, 110k sites (mac2, 5kb), 50% majority consensus tree of 100 bootstrap replicates (5M each)") +
    scale_fill_manual(values = c(`(0,50]`   = "transparent",
                                 `(50,70]`  = "white",
                                 `(70,90]`  = "gray",
                                 `(90,100]` = "black"),
                      drop = FALSE) +
    scale_size_manual(values = c(`(0,50]`   = 0,
                                 `(50,70]`  = 1.25,
                                 `(70,90]`  = 1.25,
                                 `(90,100]` = 1.25),
                      na.value = 0,
                      drop = FALSE) +
    scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
    guides(color = guide_legend(title = "Region"),
           fill = guide_legend(title = "Support", title.position = "top"),
           size = guide_legend(title = "Support", title.position = "top")
    ) +
    theme(legend.position = c(0.12, 0.6),
          legend.text = element_text(color = "gray20"),
          legend.title = element_text(color = "gray20"),
          legend.spacing.y = unit(5, "mm")
    )
)


### Save to file
ggsave(plot = t,
       filename = "figures/FigS3_svdq.pdf",
       width = 10,
       height = 14,
       device = cairo_pdf,
       bg = "transparent")



### ============================================================================
### Alternative figure versions

### Colorize only monophyletic species
# nonmono <- c("abe", "aff", "chl", "esp", "gem", "gut", "nig", "pue", "ran", "tan", "uni")
# 
# colstandard2 <- colstandard %>%
#   mutate(Color2 = case_when(Species %in% nonmono ~ "gray30", TRUE ~ Color)) %>%
#   select(-Color) %>% rename(Color = Color2)
# 
# scol <- tree %>%
#   rename(Species = spec) %>%
#   left_join(colstandard2) %>%
#   replace_na(list(Color = "gray30")) %>%   # internal branch color
#   pull(Color)
# 
# 
# ### Draw tree
# (t2 <- ggtree(tr = tree, 
#               layout = 'rectangular',
#               size = 0.5,
#               color = scol) +
#     geom_nodepoint(data = tree %>% filter(!isTip, support_class != "(0,50]"),
#                    aes(fill = support_class,
#                        size = support_class),
#                    shape = 21,
#                    color = "gray30") +
#     scale_fill_manual(values = c(`(0,50]`   = "transparent",
#                                  `(50,70]`  = "white",
#                                  `(70,90]`  = "gray",
#                                  `(90,100]` = "black"),
#                       drop = FALSE) +
#     scale_size_manual(values = c(`(0,50]`   = 0,
#                                  `(50,70]`  = 1.25,
#                                  `(70,90]`  = 1.25,
#                                  `(90,100]` = 1.25),
#                       na.value = 0,
#                       drop = FALSE) +
#     scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
#     guides(color = guide_legend(title = "Region"),
#            fill = guide_legend(title = "Support", title.position = "top"),
#            size = guide_legend(title = "Support", title.position = "top")
#     ) +
#     theme(legend.position = c(0.12, 0.6),
#           legend.text = element_text(color = "gray20"),
#           legend.title = element_text(color = "gray20"),
#           legend.spacing.y = unit(5, "mm")
#     )
# )
# 
# 
# ### Save to file
# ggsave(plot = t2,
#        filename = "figures/FigS3_svdq2.pdf",
#        width = 4,
#        height = 14,
#        device = cairo_pdf,
#        bg = "transparent")
# 
# 
# 
# ### ============================================================================
# ### Read in best SVDQ tree
# phylo <- read.nexus("results/trees/phylo2e_m2k5_svdq_1J.nex")
# 
# 
# ### Root tree
# outgroup <- c("tab", "tig", "tor")
# (outsamples <- grep(paste(outgroup, collapse = "|"), phylo$tip.label, value = TRUE))   # identify samples matching outgroup pattern
# 
# rooted <- root(phy = phylo, outgroup = outsamples, edgelabel = TRUE)
# # write.tree(rooted, file = "phylo2e_m2k5_svdq_1K.nwk")
# 
# 
# ### Add species / location labels and support categories
# (tree <- ggtree(rooted) %>%
#   .$data %>%
#   mutate(spec = if_else(isTip, str_sub(label, -6, -4), "ungrouped"),
#          loc = if_else(isTip, str_sub(label, -3, -1), "ungrouped")
#   )
# )
# 
# 
# ### Set species colors
# colstandard <- read_tsv("metadata/species_colors.tsv",
#                         col_types = "cc", 
#                         col_names = TRUE)
# 
# scol <- tree %>%
#   rename(Species = spec) %>%
#   left_join(colstandard) %>%
#   replace_na(list(Color = "gray20")) %>%   # internal branch color
#   pull(Color)
# 
# 
# ### Draw tree
# (t3 <- ggtree(tr = tree, 
#              layout = 'rectangular',
#              size = 0.5,
#              color = scol) +
#     geom_tiplab(color = "gray60", 
#                 size = 1.5, 
#                 hjust = -0.1,
#                 offset = 0.0005) +
#     # labs(title = "Hamlet phylogeny, multispecies coalescent model",
#     #      subtitle = "SVDQuartets, 110k sites (mac2, 5kb), optimal tree (490M quartets)") +
#     theme(legend.text = element_text(color = "gray20"),
#           legend.title = element_text(color = "gray20"),
#           legend.spacing.y = unit(5, "mm")
#     )
# )
# 
# 
# ### Save to file
# ggsave(plot = t3,
#        filename = "figures/FigS3_svdq3.pdf",
#        width = 10,
#        height = 14,
#        device = cairo_pdf,
#        bg = "transparent")



### ===========================================================================
### Quartetsampling support values

# phylo <- read.tree("qs/phylo2e_m2k5_svdq_1K_qs.labeled.tre.figtree") %>%
#   list_modify("edge.length" = NULL)
