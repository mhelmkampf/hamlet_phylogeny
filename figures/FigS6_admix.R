### ============================================================================
### Hamlet phylogeny: R code to plot Fig. S6 (Admixture ancestry proportions and
### GNN proportions)
### Created Nov 2023 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Read in data
prop2 <- read_delim(paste0("results/admix/AdmcvProp_phyps2e_m2n1l5_k2.tsv"), delim = " ", col_names = "Sample")
prop3 <- read_delim(paste0("results/admix/AdmcvProp_phyps2e_m2n1l5_k3.tsv"), delim = " ", col_names = "Sample")
prop4 <- read_delim(paste0("results/admix/AdmcvProp_phyps2e_m2n1l5_k4.tsv"), delim = " ", col_names = "Sample")
prop5 <- read_delim(paste0("results/admix/AdmcvProp_phyps2e_m2n1l5_k5.tsv"), delim = " ", col_names = "Sample")
prop6 <- read_delim(paste0("results/admix/AdmcvProp_phyps2e_m2n1l5_k6.tsv"), delim = " ", col_names = "Sample")
prop7 <- read_delim(paste0("results/admix/AdmcvProp_phyps2e_m2n1l5_k7.tsv"), delim = " ", col_names = "Sample")

prop2 <- prop2 %>% mutate(Species = str_sub(Sample, -6, -4), Location = str_sub(Sample, -3, -1))
prop3 <- prop3 %>% mutate(Species = str_sub(Sample, -6, -4), Location = str_sub(Sample, -3, -1))
prop4 <- prop4 %>% mutate(Species = str_sub(Sample, -6, -4), Location = str_sub(Sample, -3, -1))
prop5 <- prop5 %>% mutate(Species = str_sub(Sample, -6, -4), Location = str_sub(Sample, -3, -1))
prop6 <- prop6 %>% mutate(Species = str_sub(Sample, -6, -4), Location = str_sub(Sample, -3, -1))
prop7 <- prop7 %>% mutate(Species = str_sub(Sample, -6, -4), Location = str_sub(Sample, -3, -1))

  
### Assign species and order by region
gulf_sp <- c("atl", "cas", "eco", "esp", "flo", "gem")
caribbean_sp <- c("aff", "chl", "gum", "gut", "ind", "lib", "may", "pro", "tan")
both_sp <- c("abe", "nig", "pue", "uni", "ran")

order_sp <- c(gulf_sp, both_sp, caribbean_sp)

prop2$Species <- factor(prop2$Species, levels = order_sp)

ord_sp <- prop2[order(prop2$Species, prop2$Location), ] %>% pull(Sample)
  

### Set species colors
colstandard <- read_tsv("metadata/species_colors.tsv",
                        col_types = "cc", 
                        col_names = TRUE)
  
  
### Pivot to long format
long2 <- pivot_longer(prop2, cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion")
long3 <- pivot_longer(prop3, cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion")
long4 <- pivot_longer(prop4, cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion")
long5 <- pivot_longer(prop5, cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion")
long6 <- pivot_longer(prop6, cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion")
long7 <- pivot_longer(prop7, cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion")

  
### Plot ancestry proportions
(p2 <- ggplot(long2, aes(x = fct_relevel(Sample, ord_sp), y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = NULL,
         tag = paste0("k = 2")) +
    scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 13),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)

(p3 <- ggplot(long3, aes(x = fct_relevel(Sample, ord_sp), y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = NULL,
         tag = paste0("k = 3")) +
    scale_fill_manual(values = c("#4DAF4A", "#377EB8", "#E41A1C")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 13),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)

(p4 <- ggplot(long4, aes(x = fct_relevel(Sample, ord_sp), y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = NULL,
         tag = paste0("k = 4")) +
    scale_fill_manual(values = c("#4DAF4A", "#377EB8", "#E41A1C", "#FFFF33")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 13),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)

(p5 <- ggplot(long5, aes(x = fct_relevel(Sample, ord_sp), y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = NULL,
         tag = paste0("k = 5")) +
    scale_fill_manual(values = c("#FFFF33","#984EA3", "#4DAF4A", "#E41A1C", "#377EB8")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 13),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)

(p6 <- ggplot(long6, aes(x = fct_relevel(Sample, ord_sp), y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = NULL,
         tag = paste0("k = 6")) +
    scale_fill_manual(values = c("#FFFF33", "#E41A1C", "#4DAF4A", "#377EB8", "#984EA3", "#FF7F00")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 13),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)

(p7 <- ggplot(long7, aes(x = fct_relevel(Sample, ord_sp), y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = NULL,
         tag = paste0("k = 7")) +
    scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#FFFF33", "#377EB8", "#984EA3", "#A65628", "#FF7F00")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 13),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          # axis.text.x = element_text(hjust = 1, vjust = 1, angle = 35, size = 2,
          #                            margin = unit(c(-3, 0, 0, 0), "mm")),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
)
  
  
### Plot species bar
coord_samples <- prop2[order(prop2$Species, prop2$Location), ] %>%
  mutate(Ord_nr = row_number()) %>%
  select(Sample, Species, Ord_nr)

blocks_sp <- coord_samples %>%
  group_by(Species) %>%
  summarise(Start = min(Ord_nr) - 1,
            End = max(Ord_nr)) %>%
  left_join(colstandard) %>%
  replace(is.na(.), "gray80")

(p_sp <- ggplot(blocks_sp, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf)) +
    geom_rect(aes(fill = Species), color = "grey20", lwd = 0) +
    scale_y_continuous(breaks = .5, limits = c(0, 1)) +
    scale_x_discrete(breaks = coord_samples$Ord_nr,
                     labels = coord_samples$Species,
                     expand = c(0, 0)) +
    scale_fill_manual(values = colstandard$Color) +   ### currently wrong order
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 3, 1), "mm")
    )
)

  
### Combine admixture plots
(p_all <- plot_grid(p_sp, p2, p3, p4, p5, p6, p7, 
                    align = "v", 
                    nrow = 8,
                    rel_heights = c(0.25, 1, 1, 1, 1, 1, 1)
                    ))


### Add shared axis labels
xlabel <- textGrob("Samples",
                   y = unit(14, "mm"),
                   gp = gpar(col = "gray20", fontsize = 14, fontface = "bold"))

ylabel <- textGrob("Ancestry proportion",
                   x = unit(2, "mm"),
                   gp = gpar(col = "gray20", fontsize = 14, fontface = "bold"), rot = 90)

p_final <- grid.arrange(arrangeGrob(p_all, left = ylabel, bottom = xlabel))

  
### Save as PDF
ggsave(
  filename = "figures/FigS6_admix.pdf",
  plot = p_final,
  width = 30,
  height = 30,
  units = "cm",
  device = cairo_pdf
)



### ============================================================================
### GNN

### Read in data
# (gnn <- read_delim("results/gnn/gnn_LG02_cld.csv", delim = ",", 
#                    col_names = TRUE) %>%
#     rename(Node = `Sample node`,
#            GNN_gulf = gulf,
#            GNN_carib = carib) %>%
#     mutate(Haplotype = paste0(Individual, "-", Node),
#            Species = str_sub(Individual, -6, -4),
#            Location = str_sub(Individual, -3, -1),
#            Population = str_sub(Individual, -6, -1)) %>%
#     arrange(desc(GNN_gulf)) %>%
#     select(Haplotype, Population, Clade, GNN_gulf, GNN_carib)
# )
# 
# 
# ### Pivot to long format
# long <- pivot_longer(gnn,
#                      cols = starts_with("GNN"),
#                      names_to = "GNN_clade",
#                      values_to = "Proportion")
# 
# 
# ### Plot GNN proportions
# (g <- long %>%
#     arrange(GNN_clade, Proportion) %>%
#     mutate(Haplotype = factor(Haplotype, levels = unique(Haplotype))) %>%
#     ggplot(aes(x = Haplotype, y = Proportion, fill = GNN_clade)) +
#     geom_bar(position = "fill", stat = "identity", alpha = 0.9) +
#     labs(x = NULL,
#          y = "GNN proportions") +
#     scale_fill_manual(values = c("peru", "#837ABE")) +
#     theme_minimal(base_size = 11) +
#     theme(panel.grid.minor = element_blank(),
#           panel.grid.major.x = element_blank(),
#           plot.title = element_text(size = 13),
#           axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
#           axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
#           axis.text.x = element_text(hjust = 1, vjust = 1, angle = 35, size = 1,
#                                      margin = unit(c(-2, 0, 0, 0), "mm")),
#           # plot.tag = element_text(angle = -90),
#           # plot.tag.position = c(1.02, 0.5),
#           legend.position = "none",
#           plot.margin = unit(c(1, 10, 1, 1), "mm")
#     )
# )
# 
# 
# ### Save as PDF
# ggsave(
#   filename = "figures/FigS6_gnn.pdf",
#   plot = g,
#   width = 40,
#   height = 5,
#   units = "cm",
#   device = cairo_pdf
# )
# 
# 
# ### Combine panels (patchwork)
# # (p_comb <- p_all /
# #     g +
# #     plot_annotation(tag_levels = c("a")) +
# #     plot_layout(heights = c(5, 1)) +
# #     plot_layout(widths = c(2, 1))
# # )
# 
# 
# ### Combine plots (cowplot)
# (p_comb <- plot_grid(p_final, g,
#                      nrow = 2,
#                      rel_heights = c(0.8, 0.2)
#                      )
# )
# 
# 
# ### Save as PDF
# ggsave(
#   filename = "figures/FigS6_admixGnn.pdf",
#   plot = p_comb,
#   width = 20,
#   height = 28,
#   units = "cm",
#   device = cairo_pdf
# )



### ============================================================================
### Alternative plot versions (GNN species color bar)

# ### Create species color table
# scol <- long %>%
#   arrange(Proportion) %>%
#   mutate(spec = str_sub(Population, -6, -4),
#          position = seq(1:length(long$Haplotype)))
# 
# 
# ### Plot species colors
# (s <- ggplot(scol, aes(xmin = position, xmax = position + 1, 
#                        ymin = -Inf, ymax = Inf)) +
#     geom_rect(aes(fill = spec)) +
#     scale_x_discrete(breaks = seq(1:length(long$Haplotype))) +
#     scale_y_continuous(breaks = .5, limits = c(0, 1)) +
#     scale_fill_manual(values = colstandard$Color) +
#     labs(x = NULL, 
#          y = "Species") +
#     theme_minimal() +
#     theme(panel.grid = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           # axis.title.y = element_text(angle = 0, vjust = 0.5, 
#           #                             color = "gray20", size = 9,
#           #                             margin = margin(r = -40)),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           legend.position = "none",
#           plot.margin = unit(c(-3, 0, 3, 0), "mm")
#     )
# )
# 
# 
# ### Combine plots
# (p_comb <- plot_grid(g, s, 
#                      align = "v", 
#                      nrow = 2, 
#                      rel_heights = c(0.93, 0.07),
#                      scale = c(1, 1.015))
# )
