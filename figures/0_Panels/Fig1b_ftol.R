### ============================================================================
### Hamlet phylogeny: R code to plot Fig. 1b (Serraninae FToL phylogeny)
### Created Jan 2024 by Martin Helmkampf
### ============================================================================


### Preparations
library(ape)
library(treeio)
library(ggtree)
library(tidyverse)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Read in and root tree
phylo <- read.tree(paste("results/trees/ftol_v2.treefile"))


### Root tree
rooted <- root(phy = phylo, outgroup = "Epinephelus_maculatus", edgelabel = TRUE)


### Edit tip labels
rooted$tip.label <- rooted$tip.label %>% 
  str_replace(pattern = "PL17_88abebel", "Hypoplectrus aberrans") %>%
  str_replace(pattern = "19174affgun", "Hypoplectrus affinis") %>%
  str_replace(pattern = "54786atlliz", "Hypoplectrus atlahua") %>%
  str_replace(pattern = "54649casliz", "Hypoplectrus castroaguirrei") %>%
  str_replace(pattern = "PL17_38chlpri", "Hypoplectrus chlorurus") %>%
  str_replace(pattern = "62576ecoarc", "Hypoplectrus ecosur") %>%
  str_replace(pattern = "62558floarc", "Hypoplectrus floridae") %>%
  str_replace(pattern = "62570gemarc", "Hypoplectrus gemma") %>%
  str_replace(pattern = "23301gumboc", "Hypoplectrus gummigutta") %>%
  str_replace(pattern = "19076gutbar", "Hypoplectrus guttavarius") %>%
  str_replace(pattern = "PL17_64indpri", "Hypoplectrus indigo") %>%
  str_replace(pattern = "HypoHaiti1libhai", "Hypoplectrus liberte") %>%
  str_replace(pattern = "PL17_122maybel", "Hypoplectrus maya") %>%
  str_replace(pattern = "18906nigboc", "Hypoplectrus nigricans") %>%
  str_replace(pattern = "20650prohon", "Hypoplectrus providencianus") %>%
  str_replace(pattern = "19104puegun", "Hypoplectrus puella") %>%
  str_replace(pattern = "FL0880ransan", "Hypoplectrus randallorum") %>%
  str_replace(pattern = "PL17_72tanpri", "Hypoplectrus sp. (tan)") %>%
  str_replace(pattern = "PL17_136uniflk", "Hypoplectrus unicolor") %>%
  str_replace(pattern = "20480tabhon", "Serranus tabacarius") %>%
  str_replace(pattern = "PL17_21tigboc", "Serranus tigrinus") %>%
  str_replace(pattern = "Bocas16.3torboc", "Serranus tortugarum") %>%
#
  str_replace(pattern = "_", " ")
  # str_replace(pattern = "Za.", "Pl.") %>%
  # str_replace(pattern = "Hy.", "H.") %>% 
  # str_replace(pattern = "Di.", "D.") %>% 
  # str_replace(pattern = "Ep.", "E.")


### Identify long branches (requires ggtree object index)
tree_data <- ggtree(rooted)$data
longg <- which(tree_data$branch.length > mean(tree_data$branch.length) * 10)


### Compress long branches (requires rooted phylo object index)
# longr <- which(rooted$edge.length > mean(rooted$edge.length) * 10)
# rooted$edge.length[longr] <- rooted$edge.length[longr] * 0.5


### Define groups
hamlets <- rooted$tip.label[grepl(pattern = "Hypoplectrus", rooted$tip.label)]
our_taxa <- c(hamlets, "Serranus tabacarius", "Serranus tigrinus", "Serranus tortugarum")


### Add species / location labels and support categories
(tree <- ggtree(rooted, 
                layout = "rectangular", 
                ladderize = TRUE, 
                right = TRUE) %>%
    .$data %>%
    mutate(support = as.numeric(if_else(!isTip, label, "NA")),
           support_class = cut(as.numeric(support), c(0,50,70,90,100)) %>%
             as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"))) %>%
           # branch_type = case_when(node %in% longg ~ "broken", 
           #                         TRUE ~ "whole")) %>%
    groupOTU(.node = hamlets)
)


### Define neutral color
clr_neutral <- rgb(0.2, 0.2, 0.2)


### Draw tree
(t <- ggtree(tree,
             aes(color = clr_neutral,
                 # linetype = branch_type
                 )) +
  geom_tiplab(aes(color = group, label = if_else(label %in% our_taxa, str_c(label,"*"), label)),   # add asterisks to our taxa
              size = 3.5, 
              hjust = -0.05,
              fontface = "italic") +
  ggplot2::xlim(0, 0.19) +   # add extra space for long labels
  geom_nodepoint(aes(fill = support_class,
                     size = support_class),
                     shape = 21) +
  # labs(title = "Serraninae, updated FToL phylogeny",
  #     subtitle = "23 nuclear and mt genes, partition model, IQ-TREE, 200 non-parametric bootstraps") +
  scale_color_manual(values = c(clr_neutral, clr_neutral, "firebrick")) +
  scale_fill_manual(values = c(`(0,50]`   = "transparent",
                               `(50,70]`  = "white",
                               `(70,90]`  = "gray",
                               `(90,100]` = "black"),
                    breaks = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"),
                    drop = FALSE) +
  scale_size_manual(values = c(`(0,50]`   = 0,
                               `(50,70]`  = 1.5,
                               `(70,90]`  = 1.5,
                               `(90,100]` = 1.5),
                    breaks = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"),
                    na.value = 0,
                    drop = FALSE) +
  scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
  geom_treescale(color = "gray", 
                 fontsize = 3, 
                 x = 0, 
                 y = 15) +
  guides(color = "none",
         fill = guide_legend(title = "Support", title.position = "top"),
         size = guide_legend(title = "Support", title.position = "top")
         ) +
  theme_void() +
  theme(plot.margin =  margin(l = 10),
        plot.title = element_text(size = 13, face = "bold", vjust = 3, margin = margin(15, 0, 5, 0)),
        plot.subtitle = element_text(size = 11, margin = margin(0, 0, 30, 0)),
        legend.position = c(0.05,0.05),
        legend.justification = c(0,0),
        legend.title.align = 0,
        legend.text = element_text(color = clr_neutral),
        legend.title = element_text(color = clr_neutral))
  )


### Save to file
ggsave(plot = t,
       filename = "figures/Fig1b_ftol.pdf",
       width = 8,
       height = 6.5,
       device = cairo_pdf,
       bg = "transparent")
