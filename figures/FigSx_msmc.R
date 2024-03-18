### ============================================================================
### Hamlet phylogeny: R code to plot Fig. Sx (Ne and cross coalescence rates)
### Created Mar 2024 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)
library(patchwork)
library(GenomicOriginsScripts)   # see https://k-hench.github.io/GenomicOriginsScripts/index.html

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list=ls())
dev.off()


### Set paths
Ne_path <- "results/msmc/Ne/"
cc_path <- "results/msmc/cc/"


### Read in grouping files
Ne_grouping <- read_tsv("data/msmc/msmc_grouping_phylo2e-n3.tsv")
cc_grouping <- read_tsv("data/msmc/cc_phylo2e_gom.tsv")


### Locate results
Ne_files <- dir(Ne_path, pattern = "*.final.txt")
cc_files <- dir(cc_path, pattern = "*.final.txt")


### Import results
Ne_data <- Ne_files %>%
  map_dfr(.f = get_msmc, msmc_path = Ne_path)

cc_data <- cc_files %>%
  map_dfr(get_cc, cc_groups = cc_grouping, cc_path = cc_path)


### Set species colors
colstandard <- read_tsv("metadata/species_colors.tsv", 
                        col_types = "cc", 
                        col_names = TRUE) %>%
  filter(Species != "cas")


### Define groups
small <- c("atl", "eco", "flo")
large <- c(subset(colstandard$Species, !(colstandard$Species %in% small)), "cas")


### Plot Ne
(p_Ne <- Ne_data %>%
  filter(!time_index %in% c(0:4, 27:31)) %>%   # remove the two first and last time segments
  ggplot(aes(x = YBP, y = Ne, group = run_nr, color = spec)) +
  geom_line(linewidth = 0.5) +
  annotation_logticks(sides = "tl", color = "lightgray", size = plot_lwd) +   # add guides for the logarithmic axes
  scale_color_manual(values = colstandard$Color, label = colstandard$Species) +
  scale_x_log10(expand = c(0, 0),
                breaks = c(10^3, 10^4, 10^5),
                position = "top",
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = c(10^3, 10^4, 10^5, 10^6),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Generations Before Present",
       y = expression(Effective~Population~Size~(italic(N[e])))) +
  coord_cartesian(xlim = c(500, 3*10^5)) +
  guides(color = guide_legend(title = "Species", 
                              title.position = "top",
                              byrow = TRUE)) +
  theme_minimal(base_size = 12) +
  theme(axis.ticks = element_line(colour = "lightgray"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.spacing.y = unit(10, "pt")
  )
)


### Filter cross-coalescence rates
cc_filt <- cc_data %>%
  filter(!time_index %in% c(0:4, 27:31)) %>%   # remove the two first and last time segments
  mutate(comp = case_when(spec_1 %in% small & spec_2 %in% small ~ "small",
                          spec_1 %in% large & spec_2 %in% large ~ "large",
                          TRUE ~ "between"))
 

### Plot cross-coalescence rates
(p_cc <- cc_filt %>%
  ggplot(aes(x = YBP, y = Cross_coal, group = run_nr, color = comp)) +
  geom_line(linewidth = 0.6, alpha = 0.3) +
  annotation_logticks(sides = "b", color = "lightgray", size = plot_lwd) +
  scale_color_manual(values = c("darkgreen", "peru", "purple")) +
  scale_x_log10(expand = c(0, 0),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = NULL,
       y = "Cross-coalescence Rate") +
  coord_cartesian(xlim = c(500, 3*10^5)) +
  guides(color = guide_legend(title = "Comparison", 
                              title.position = "top",
                              byrow = TRUE)) +
  theme_minimal(base_size = 12) +
  theme(axis.ticks = element_line(colour = "lightgray"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.spacing.y = unit(10, "pt")
    )
)


### Combine panels
(p_comb <- p_Ne /
  p_cc  +
  plot_annotation(tag_levels = c("a")) &
  theme(panel.grid.major = element_line(linewidth = plot_lwd)
  )
)


### Save plot to file
ggsave(plot = p_comb,
       filename = "figures/FigSx_msmc.pdf",
       width = 7,
       height = 10,
       device = cairo_pdf,
       bg = "transparent")



### ============================================================================
### Alternative figure versions

### Between clades only, by small clade
# cc_alt1 <- cc_data %>%
#   filter(!time_index %in% c(0:4, 27:31)) %>%   # remove the two first and last time segments
#   mutate(comp = case_when(spec_1 %in% small & spec_2 %in% small ~ "small",
#                           spec_1 %in% large & spec_2 %in% large ~ "large",
#                           TRUE ~ "between")) %>%
#   filter(comp == "between") %>%
#   mutate(cont = case_when(spec_1 == "atl" | spec_2 == "atl" ~ "includes atl",
#                           spec_1 == "eco" | spec_2 == "eco" ~ "includes eco",
#                           spec_1 == "flo" | spec_2 == "flo" ~ "includes flo",
#                           TRUE ~ "NA")
#          )
# 
# 
# (p_alt1 <- cc_alt1 %>%
#     ggplot(aes(x = YBP, y = Cross_coal, group = run_nr, color = cont)) +
#     geom_line(linewidth = 0.6, alpha = 0.5) +
#     annotation_logticks(sides = "b", color = "lightgray", size = plot_lwd) +
#     scale_x_log10(expand = c(0, 0),
#                   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#     labs(x = NULL,
#          y = "Cross-coalescence Rate") +
#     coord_cartesian(xlim = c(500, 3*10^5)) +
#     guides(color = guide_legend(title = "Between only", 
#                                 title.position = "top",
#                                 byrow = TRUE)) +
#     theme_minimal(base_size = 12) +
#     theme(axis.ticks = element_line(colour = "lightgray"),
#           panel.grid.minor.x = element_blank(),
#           panel.grid.minor.y = element_blank(),
#           legend.spacing.y = unit(10, "pt")
#     )
# )
# 
# 
# ggsave(plot = p_alt1,
#        filename = "figures/FigSx_msmc_alt1.pdf",
#        width = 7,
#        height = 5,
#        device = cairo_pdf,
#        bg = "transparent")
# 
# 
# ### Small clade only, by comparison
# cc_alt2 <- cc_data %>%
#   filter(!time_index %in% c(0:4, 27:31)) %>%   # remove the two first and last time segments
#   mutate(comp = case_when(spec_1 %in% small & spec_2 %in% small ~ "small",
#                           spec_1 %in% large & spec_2 %in% large ~ "large",
#                           TRUE ~ "between")) %>%
#   filter(comp == "small")
# 
# 
# (p_alt2 <- cc_alt2 %>%
#     ggplot(aes(x = YBP, y = Cross_coal, group = run_nr, color = run)) +
#     geom_line(linewidth = 0.6, alpha = 0.5) +
#     annotation_logticks(sides = "b", color = "lightgray", size = plot_lwd) +
#     scale_color_manual(labels = c("atl -- eco", "atl -- flo", "eco -- flo"),
#                        values = c("red", "green", "blue")) +
#     scale_x_log10(expand = c(0, 0),
#                   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#     labs(x = NULL,
#          y = "Cross-coalescence Rate") +
#     coord_cartesian(xlim = c(500, 3*10^5)) +
#     guides(color = guide_legend(title = "Small only", 
#                                 title.position = "top",
#                                 byrow = TRUE)) +
#     theme_minimal(base_size = 12) +
#     theme(axis.ticks = element_line(colour = "lightgray"),
#           panel.grid.minor.x = element_blank(),
#           panel.grid.minor.y = element_blank(),
#           legend.spacing.y = unit(10, "pt")
#     )
# )
# 
# 
# ### Save plot to file
# ggsave(plot = p_alt2,
#        filename = "figures/FigSx_msmc_alt2.pdf",
#        width = 7,
#        height = 5,
#        device = cairo_pdf,
#        bg = "transparent")
