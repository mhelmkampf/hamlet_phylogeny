### ============================================================================
### Hamlet phylogeny: R code to plot Fig. 2a, concatenated SNPs tree
### Created 24 Oct 2023 by Martin Helmkampf
### Modified by Floriane Coulmance for final figures
### usage:
### Rscript Fig2a_iqtree.R
### ============================================================================


### Clear workspace
dev.off()
rm(list=ls())


### Load needed library
library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(tidyverse)
library(ggnewscale)
library(ggtreeExtra)
library(ggtext)


### Read in tree, change to your path
phylo <- read.tree("/hamlet_phylogeny/results/trees/iqtree_phylo2e_m2k5_GTRL_1A.treefile")



# 1) Data preparation for plotting
# -------------------------------------------------------------------------------------------------------------------

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


### Define groups
gulf <- c("liz", "tam", "ala", "arc", "are", "flk")
caribbean <- c("bel", "boc", "gun", "hon", "san", "qui")
atlantic <- c("bar", "hai", "pri")



# 2) Graphic preparation for plotting
# -------------------------------------------------------------------------------------------------------------------

### Add species / location labels and support categories
(tree <- ggtree(rooted) %>%
  .$data %>%
  mutate(spec = if_else(isTip, str_sub(label, -6, -4), "ungrouped"),
         loc = if_else(isTip, str_sub(label, -3, -1), "ungrouped"),
         region = case_when(
           loc %in% gulf ~ "Gulf of Mexico",
           loc %in% caribbean ~ "Caribbean",
           loc %in% atlantic ~ "Atlantic (incl. Florida)",
           # loc %in% florida ~ "Florida",
           TRUE ~ "NA"),
         support = as.numeric(if_else(!isTip, label, "NA")),
         support_class = cut(as.numeric(support), c(0, 50, 70, 90, 100)) %>%
           as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")),
         branch_type = case_when(node %in% longg ~ "broken", 
                                 TRUE ~ "whole")
  )
)


### Set species colors
colstandard <- read_tsv("/hamlet_phylogeny/metadata/species_colors.tsv", 
                        col_types = "cc", 
                        col_names = TRUE)

scol <- tree %>%
  rename(Species = spec) %>%
  left_join(colstandard) %>%
  replace_na(list(Color = "gray20")) %>%   # internal branch color
  pull(Color)


### Set species legend
logos_spec <- c(abe = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_aberrans.l.cairo.png' width='80' /><br>*H. aberrans*",
                chl = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_chlorurus.l.cairo.png' width='80' /><br>*H. chlorurus*",
                flo = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_floridae.l.cairo.png' width='80' /><br>*H. floridae*",
                gem = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_gemma.l.cairo.png' width='80' /><br>*H. gemma*",
                gum = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_gumigutta.l.cairo.png' width='80' /><br>*H. gummigutta*",
                gut = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_guttavarius.l.cairo.png' width='80' /><br>*H. guttavarius*",
                ind = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_indigo.l.cairo.png' width='80' /><br>*H. indigo*",
                may = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_maya.l.cairo.png' width='80' /><br>*H. maya*",
                nig = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_nigricans.l.cairo.png' width='80' /><br>*H. nigricans*",
                pue = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_puella.l.cairo.png' width='80' /><br>*H. puella*",
                ran = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_randallorum.l.cairo.png' width='80' /><br>*H. randallorum*",
                aff = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_affinis.l.cairo.png' width='80' /><br>*H. affinis*",
                uni = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_unicolor.l.cairo.png' width='80' /><br>*H. unicolor*",
                tan = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_sp1.l.cairo.png' width='80' /><br>*H. sp1*",
                pro = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_providencianus.l.cairo.png' width='80' /><br>*H. providencianus*",
                atl = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_atlahua.l.cairo.png' width='80' /><br>*H. atlahua*",
                cas = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_castroaguirrei.l.cairo.png' width='80' /><br>*H. castroaguirrei*",
                eco = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_ecosur.l.cairo.png' width='80' /><br>*H. ecosur*",
                lib = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_liberte.l.cairo.png' width='80' /><br>*H. liberte*",
                tor = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/S_tortugarum.l.cairo.png' width='80' /><br>*S. tortugarum*",
                tab = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/S_tabacarius.l.cairo.png' width='80' /><br>*S. tabacarius*",
                tig = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/S_tigrinus.l.cairo.png' width='80' /><br>*S. tigrinus*",
                esp = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_espinoza.l.cairo.png' width='80' /><br>*H. sp2*",
                ungrouped = "*not a clade*")



# 3) Plotting
# -------------------------------------------------------------------------------------------------------------------

### Plot tree
(t <- ggtree(tr = tree, 
             layout = 'fan',
             size = 1.2,
             color = scol,
             aes(linetype = branch_type)) +
    geom_tippoint(aes(color = spec),
                  size = 3,
                  alpha = 0.5) +
    geom_nodepoint(data = tree %>% filter(!isTip, support_class != "(0,50]"),
                   aes(fill = support_class,
                       size = support_class),
                   shape = 21,
                   color = "gray20") +
    scale_color_manual(values = c(
      abe = "#DFDF8D",
      gum = "#F99729",
      nig = "#333333",
      pue = "#E48175",
      aff = "#FFB1D8",
      uni = "#B3B3B3",
      flo = "#A155E7",
      chl = "#8B4513",
      gem = "#1E90FF",
      gut = "#FFEA00",
      ind = "#22198E",
      may = "#7196AE",
      ran = "#8AC3BA",
      tan = "#D2B48C",
      atl = "#5C7A1E",
      cas = "#E3280D",
      eco = "#FFBF70",
      lib = "#C3A8DC",
      pro = "#8BCF4C",
      esp = "#E93E95"),
    labels = logos_spec,
    guide = "none") +
  scale_fill_manual(values = c(`(0,50]`   = "transparent",
                                 `(50,70]`  = "white",
                                 `(70,90]`  = "gray",
                                 `(90,100]` = "black"),
                      drop = FALSE,
                    guide = "none") +
    scale_size_manual(values = c(`(0,50]`   = 0,
                                 `(50,70]`  = 4,
                                 `(70,90]`  = 4,
                                 `(90,100]` = 6),
                      na.value = 0,
                      drop = FALSE) +
    scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
    guides(color = guide_legend(title = "Species", nrow = 5, title.position = "top", order = 1),
           fill = guide_legend(title = "Support", title.position = "top", order = 2),
           size = guide_legend(title = "Support", title.position = "top", order = 2)
           ) +
    ggtree::geom_treescale(width = 0.01,
                           # offset = -4, 
                           fontsize = 3, 
                           color = "gray60",
                           x = 0.05, y = 0) +
    theme(legend.position = "bottom", #c(0.12, 0.3),
          plot.title = element_text(size = 14, color = "gray20", face = "bold", vjust = 3, margin = margin(t = 10)),
          plot.subtitle = element_text(size = 10, color = "gray20", margin = margin(t = 5, b = 15)),
          legend.text = element_markdown(color = "gray20", size=11),
          legend.title = element_text(color = "gray20"),
          legend.spacing.y = unit(5, "mm"),
          legend.background = element_rect(fill='transparent'))
)



### Set up the outter ring heatmap colors for species and geographic site
v <- tree %>% filter(isTip == TRUE)
colrspec <- data.frame("spec" = v[,c("spec")])
rownames(colrspec) <- tree$label[tree$isTip == TRUE]

r <- v %>%
  rename(Species = spec) %>%
  left_join(colstandard) %>%
  replace_na(list(Color = "gray20"))

r <- r %>%
  mutate(
    noham=case_when(
      endsWith(region, "Caribbean") ~ "group1",
      endsWith(region, "Atlantic (incl. Florida)") ~ "group1",
      endsWith(region, "Gulf of Mexico") ~ "group1"),
    noham2=case_when(
      endsWith(region, "Caribbean") ~ "Caribbean",
      endsWith(region, "Atlantic (incl. Florida)") ~ "Atlantic",
      endsWith(region, "Gulf of Mexico") ~ "GulfofMexico")) %>%
  reframe(label, Species, Color, region, isTip, noham, noham2) %>%
  setNames(c("Label", "hamlets", "color", "location_region", "tip", "noham", "noham2"))


### Plot the final tree with outter ring heatmap for biogeographic region 
(p <- t + new_scale_fill() +
    geom_fruit(data=r, geom=geom_tile,
               mapping=aes(y=Label, x=noham, fill=noham2), alpha = 0.8, offset = 0.02, width = 0.003) + 
    scale_fill_manual(values = c(GulfofMexico = "coral2",
                                 Caribbean = "royalblue1",
                                 Florida = "olivedrab4",
                                 Atlantic = "olivedrab4"),
                      guide = "none") 
)


### Save pdf figure
ggsave(plot = p,
       filename = "/hamlet_phylogeny/figures/Fig2a_iqtree.pdf",
       width = 32,
       height = 36,
       device = cairo_pdf,
       bg = "transparent",
       limitsize = FALSE)