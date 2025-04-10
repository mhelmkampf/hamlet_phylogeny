#!/usr/bin/env Rscript
# by: Floriane Coulmance: 19/05/2023
# usage:
# Rscript Fig3c_S10-19.R <base_dir> <file_path>
# -------------------------------------------------------------------------------------------------------------------
# <base_dir> in: $BASE_DIR aka /hamlet_phylogeny/
# <file_path> path to any *.raxml.support file created by pipelines
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed libraries
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(ape)
library(ggtree)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(tidygraph)
library(ggraph)
library(patchwork)
library(ggbreak)
library(treeio)
library(phangorn)
library(ggtreeExtra)
library(ggnewscale)



# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

base_dir <- as.character(args[1]) # path to the root directory
# base_dir <- "/Users/fco/Desktop/PhD/2_CHAPTER2/chapter2/"

phyps2 <- as.character(args[2]) # path to the tree file of interest
# phyps2 <- "/Users/fco/Desktop/PhD/2_CHAPTER2/chapter2/figures/gxp_clades/large/8_regions.raxml.support"



# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

# no functions

# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


  # Read the tree file data in table
raxml_tree <- read.tree(phyps2)


# 1) Data preparation for plotting
# -------------------------------------------------------------------------------------------------------------------

  # Change incorrect labels
raxml_tree$tip.label <- raxml_tree$tip.label %>% 
  str_replace(pattern = "54761nigliz", "54761atlliz") %>%
  str_replace(pattern = "52989nigtam", "52989atltam") %>%
  str_replace(pattern = "27936flotam", "27936atltam") %>%
  str_replace(pattern = "PL17_35puepri", "PL17_35indpri") %>%
  str_replace(pattern = "19294tangun", "19294affgun") %>%
  str_replace(pattern = "62585uniarc", "62585esparc") %>%
  str_replace(pattern = "62571uniarc", "62571esparc") %>%
  str_replace(pattern = "62555uniarc", "62555esparc") %>%
  str_replace(pattern = "27698gutqui", "27698abequi")

  # Root the tree with outgroup species 
outgroup <- c("tig", "tab", "tor")
outsamples <- grep(paste(outgroup, collapse = "|"), raxml_tree$tip.label, value = TRUE)
rooted <- root(phy = raxml_tree, outgroup = outsamples, edgelabel = TRUE)

  # In the 2 next steps, here are the settings used to identify long branches and reduce them
  # Each tree was different and for the sake of visualization, those parameters could not be constant
  #                     branch_length >   |   reduced by  
  # Fig3c_LG12_4_tree:  median*55         |   5% (0.05) 
  # FigS10_casz1_wg:    median*35         |   5% (0.05)
  # FigS11_casz1_cds:   median*655        |   5% (0.05)
  # FigS12_LG04_1:      median*15         |   5% (0.05)
  # FigS13_LG08_2:      median*75         |   2% (0.02)
  # FigS14_LG08_3:      median*55         |   5% (0.05)
  # FigS15_LG12_5:      median*55         |   5% (0.05)
  # FigS16_LG17_6:      median*35         |   5% (0.05)
  # FigS17_LG19_7:      median*35         |   5% (0.05)
  # FigS18_LG23_8:      median*35         |   5% (0.05)
  # FigS19_8_reg:       median*55         |   0.5% (0.005)
  # Identify long branches
tree_data <- ggtree(rooted)$data
longg <- which(tree_data$branch.length > median(tree_data$branch.length) * 55)

  # Compress long branches
longr <- which(rooted$edge.length > median(rooted$edge.length) * 55)
rooted$edge.length[longr] <- rooted$edge.length[longr] * 0.005


# 2) Graphic preparation for plotting
# -------------------------------------------------------------------------------------------------------------------

  # Set up markdown legend
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

  # Set the different location groups
gulf <- c("liz", "tam", "ala", "arc", "are")
caribbean <- c("bel", "boc", "gun", "hon", "san", "qui")
atlantic <- c("bar", "hai", "pri", "flk")

  # Add columns to be used for graphics in the tree data 
tree <- ggtree(rooted) %>%
  .$data %>%
  mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
         geo = ifelse(isTip, str_sub(label, -3, -1), "ungrouped"),
         support = as.numeric(if_else(!isTip, label, "NA")),
         support_class = cut(support, c(0,50,70,90,100)) %>%
           as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")),
         species = case_when(
           endsWith(spec, "nig") ~ "nigricans",
           endsWith(spec, "uni") ~ "unicolor",
           endsWith(spec, "pue") ~ "puella",
           endsWith(spec, "abe") ~ "aberrans",
           endsWith(spec, "gum") ~ "gummiguta",
           endsWith(spec, "aff") ~ "affinis",
           endsWith(spec, "flo") ~ "floridae",
           endsWith(spec, "tan") ~ "tan",
           endsWith(spec, "ran") ~ "randallorum",
           endsWith(spec, "chl") ~ "chlorurus",
           endsWith(spec, "gem") ~ "gemma",
           endsWith(spec, "gut") ~ "guttavarius",
           endsWith(spec, "ind") ~ "indigo",
           endsWith(spec, "may") ~ "maya",
           endsWith(spec, "atl") ~ "atlahua",
           endsWith(spec, "cas") ~ "castroaguirrei",
           endsWith(spec, "eco") ~ "ecosur",
           endsWith(spec, "lib") ~ "liberte",
           endsWith(spec, "pro") ~ "providencianus",
           endsWith(spec, "tig") ~ "S. tigrinus",
           endsWith(spec, "tab") ~ "S. tabacarius",
           endsWith(spec, "tor") ~ "S. tortugarum",
           endsWith(spec, "esp") ~ "espinozaperezei"),
         region = case_when(geo %in% gulf ~ "GulfofMexico",
                            geo %in% caribbean ~ "Caribbean",
                            geo %in% atlantic ~ "Atlantic",
                            TRUE ~ "NA"),
         branch_type = case_when(node %in% longg ~ "broken",
                                 TRUE ~ "whole")
         )

  # Read species color codes
colstandard <- read_tsv(paste0(base_dir,"/metadata/species_colors.tsv"), 
                        col_types = "cc", 
                        col_names = TRUE)

  # Incorporate species color code in tree data
scol <- tree %>%
  mutate(Species = spec) %>%
  left_join(colstandard) %>%
  replace_na(list(Color = "gray20")) %>%   # internal branch color
  pull(Color)


# 3) Plotting
# -------------------------------------------------------------------------------------------------------------------

  # 1st tree
(p_tree <- ggtree(tr = tree,
                 size = 1.2,
                 layout = "fan",
                 color = scol,
                 aes(linetype = branch_type)) +
  geom_tippoint(aes(color=spec),
                  size = 4,
                  alpha = 0.5) + 
  geom_nodepoint(data = tree %>% filter(!isTip, support_class != "(0,50]"),
                 aes(fill = support_class,
                     size = support_class),
                 shape = 21,
                 color = "gray20") +
  scale_color_manual(values = c(#ungrouped = rgb(.6, .6, .6),
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
                                esp = "#E93E95"
                                # GulfofMexico = "coral2",
                                # Caribbean = "royalblue1",
                                # Atlantic = "olivedrab4",
                                # Florida = "olivedrab4"
                                  ),
                     labels = logos_spec,
                     guide = "none") +
  scale_fill_manual(values = c(`(0,50]` = "transparent",
                               `(50,70]` = "white",
                               `(70,90]` = "gray",
                               `(90,100]` = "black"),
                    drop=FALSE,
                    guide = "none") +
  scale_size_manual(values = c(`(0,50]` = 0,
                               `(50,70]` = 4,
                               `(70,90]` = 4,
                               `(90,100]` = 6),
                    na.value = 0,
                    drop = FALSE,
                    guide = "none") +
  scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
  ggtree::geom_treescale(width = 0.05,
                         fontsize = 3,
                         color = "gray60",
                         x = .00003, y = 0) +
  theme(legend.position = "bottom", #c(0.3, 0.5),
        plot.title = element_text(size = 14, color = "gray20", face = "bold", vjust = 3, margin = margin(t = 10)),
        plot.subtitle = element_text(size = 10, color = "gray20", margin = margin(t = 5, b = 15)),
        legend.text = element_markdown(color = "gray20", size = 11),
        legend.title = element_text(color = "gray20"),
        legend.spacing.y = unit(5, "mm"),
        legend.background = element_rect(fill='transparent'))
)

  # Set up the outter ring heatmap colors for species and geographic site
v <- tree %>% filter(isTip == TRUE)
colrspec <- data.frame("spec" = v[,c("spec")])
rownames(colrspec) <- tree$label[tree$isTip == TRUE]

r <- v %>%
  mutate(Species = spec) %>%
  left_join(colstandard) %>%
  replace_na(list(Color = "gray20"))

r <- r %>%
  mutate(
         noham=case_when(
           endsWith(Species, "nig") ~ "group1",
           endsWith(Species, "uni") ~ "group1",
           endsWith(Species, "pue") ~ "group1",
           endsWith(Species, "abe") ~ "group1",
           endsWith(Species, "gum") ~ "group1",
           endsWith(Species, "aff") ~ "group1",
           endsWith(Species, "flo") ~ "group1",
           endsWith(Species, "tan") ~ "group1",
           endsWith(Species, "ran") ~ "group1",
           endsWith(Species, "chl") ~ "group1",
           endsWith(Species, "gem") ~ "group1",
           endsWith(Species, "gut") ~ "group1",
           endsWith(Species, "ind") ~ "group1",
           endsWith(Species, "may") ~ "group1",
           endsWith(Species, "atl") ~ "group1",
           endsWith(Species, "cas") ~ "group1",
           endsWith(Species, "eco") ~ "group1",
           endsWith(Species, "lib") ~ "group1",
           endsWith(Species, "pro") ~ "group1",
           endsWith(Species, "tig") ~ "group3",
           endsWith(Species, "tab") ~ "group3",
           endsWith(Species, "tor") ~ "group3",
           endsWith(Species, "esp") ~ "group1")) %>%
  reframe(label, Species, Color, region, isTip, noham) %>%
  setNames(c("Label", "hamlets", "color", "location_region", "tip", "noham"))

  # plot the final tree with outter ring heatmap for species 
(p <- p_tree + new_scale_fill() +
  geom_fruit(data=r, geom=geom_tile,
             mapping=aes(y=Label, x=noham, fill=hamlets), offset = .03, pwidth = 0.1) + 
  scale_fill_manual(values = c(
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
               esp = "#E93E95",
               tig = "white",
               tor = "white",
               tab = "white"),
                 # GulfofMexico = "coral2",
               # Caribbean = "royalblue1",
               # Atlantic = "olivedrab4",
               # Florida = "olivedrab4"
             guide = "none") 
)

  # Save the tree figure
hypo_save(plot = p,
          filename = paste0(base_dir,"/figures/gxp_clades/large/8_regions.pdf"),
          width = 46,
          height = 50,
          device = cairo_pdf,
          bg = "transparent",
          limitsize = FALSE)

