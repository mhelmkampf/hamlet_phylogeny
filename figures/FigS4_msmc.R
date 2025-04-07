### ============================================================================
### Hamlet phylogeny: R code to plot Fig. S4 (Ne and cross coalescence rates)
### Created Mar 2024 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)
library(patchwork)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Set paths
Ne_path <- "results/msmc/Ne/"
cc_path <- "results/msmc/cc/"


### Read in grouping files
Ne_grouping <- read_tsv("data/msmc/Ne_grouping_phylo2e_n3.tsv")
cc_grouping <- read_tsv("data/msmc/cc_grouping_phylo2e_gom.tsv")


### Locate results
Ne_files <- dir(Ne_path, pattern = "*.final.txt")
cc_files <- dir(cc_path, pattern = "*.final.txt")


### Functions to convert results (from GenomicOriginsScripts package by Kosmas Hench)
get_msmc <- function (file, msmc_path, mu = 3.7e-08, gen = 1) 
{
  vroom::vroom(stringr::str_c(msmc_path, file), delim = "\t") %>% 
    tidyr::gather("Side", "time_value", 2:3) %>% dplyr::arrange(time_index) %>% 
    dplyr::mutate(YBP = time_value/mu * gen, Ne = (1/lambda)/mu, 
                  run_nr = stringr::str_replace(file, pattern = "run([0-9]*).*", 
                                                replacement = "\\1") %>% as.numeric(), spec = stringr::str_replace(file, 
                                                                                                                   pattern = "run[0-9]*\\.([a-z]*).*", replacement = "\\1"), 
                  loc = stringr::str_replace(file, pattern = "run[0-9]*\\.[a-z]{3}\\.([a-z]*).*", 
                                             replacement = "\\1"), run = stringr::str_c(spec, 
                                                                                        loc))
}

get_cc <- function (file, cc_groups, cc_path, mu = 3.7e-08, gen = 1) 
{
  cc_run <- stringr::str_replace(file, pattern = "cc_run\\.([0-9]*).*", 
                                 replacement = "\\1") %>% as.numeric()
  specs <- c(cc_groups$spec_1[cc_groups$run_nr == cc_run], 
             cc_groups$spec_2[cc_groups$run_nr == cc_run]) %>% sort()
  loc <- cc_groups$geo[cc_groups$run_nr == cc_run]
  vroom::vroom(stringr::str_c(cc_path, file), delim = "\t") %>% 
    tidyr::gather("Side", "time_value", 2:3) %>% dplyr::arrange(time_index) %>% 
    dplyr::mutate(YBP = time_value/mu * gen, Cross_coal = 2 * 
                    lambda_01/(lambda_00 + lambda_11), run_nr = cc_run, 
                  spec_1 = specs[1], spec_2 = specs[2], loc = loc, 
                  run = str_c(spec_1, loc, "-", spec_2, loc))
}


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


# Define replacement labels
spn_labels <- colstandard$Species
spn_labels[which(spn_labels == "tan")] <- "sp1"


### Define groups
small <- c("atl", "eco", "flo")
large <- c(subset(colstandard$Species, !(colstandard$Species %in% small)), "cas")


### Plot Ne
(p_Ne <- Ne_data %>%
  filter(!time_index %in% c(0:4, 27:31)) %>%   # remove the first five and last five time segments
  ggplot(aes(x = YBP, y = Ne, group = run_nr, color = spec)) +
  geom_line(linewidth = 0.5) +
  annotation_logticks(sides = "tl", color = "lightgray", linewidth = 0.5) +   # add guides for the logarithmic axes
  scale_color_manual(values = colstandard$Color, labels = spn_labels) +
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
  annotation_logticks(sides = "b", color = "lightgray", linewidth = 0.5) +
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
  plot_annotation(tag_levels = c("A")) &
  theme(panel.grid.major = element_line(linewidth = 0.5),
        plot.tag = element_text(face = "bold")
  )
)


### Save plot to file
ggsave(plot = p_comb,
       filename = "figures/FigS4_msmc.pdf",
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
