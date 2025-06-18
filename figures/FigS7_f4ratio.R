### ============================================================================
### Hamlet phylogeny: R code to plot f4 ratios and Patterson's D (Dstats trios)
### Created Nov 2023 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)
library(gtools)
library(patchwork)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Read in D-statistics
(dstats_BBAA <- read_delim("results/dstats/phylo2e_m2n1l5.dtrios_BBAA.txt", delim = "\t") %>%
  rename(Z_score = `Z-score`, p_value = `p-value`, f4_ratio = `f4-ratio`) %>%
  mutate(correction = p.adjust(p_value, method = "bonferroni")) %>%
  arrange(correction) %>%
  rename(p_adjusted = `correction`)
)

(dstats_Dmin <- read_delim("results/dstats/phylo2e_m2n1l5.dtrios_Dmin.txt", delim = "\t") %>%
    rename(Z_score = `Z-score`, p_value = `p-value`, f4_ratio = `f4-ratio`) %>%
    mutate(correction = p.adjust(p_value, method = "bonferroni")) %>%
    arrange(correction) %>%
    rename(p_adjusted = `correction`)
)


### Plot distribution of D, f4 ratio, and p
plot(dstats_BBAA$Dstatistic, ylab = "D", xlab = "trio number")
plot(dstats_BBAA$f4_ratio, ylab = "f4 ratio", xlab = "trio number")
plot(dstats_BBAA$p_adjusted, ylab = "adjusted p", xlab = "trio number")


### Declare populations
samples <- read_tsv("metadata/ids_phylo2e.txt", col_names = "Sample") %>%
  mutate(Population = str_sub(Sample, -6, -1)) %>%
  filter(!grepl("tor|tab|tig", Population))

populations <- levels(factor(samples$Population))


### Order populations
alpha <- sort(populations)

pop_by_region <- c(
  populations[str_detect(populations, "liz|tam")],
  populations[str_detect(populations, "ala|arc|are|flk")],
  populations[str_detect(populations, "bar|hai|pri")],
  populations[str_detect(populations, "bel|boc|gun|hon|san|qui")]
)


### Function to add pair label
f_sorted_pair <- function(P2, P3, ...){
  n1 <- as.numeric(factor(P2, levels = alpha))
  n2 <- as.numeric(factor(P3, levels = alpha))
  char_sorted <- if (n1 < n2) {c(P2, P3)} else {c(P3, P2)}
  str_c(char_sorted[[1]], "-", char_sorted[[2]])
}


### Find max/min values for D, p, f4
(minmax_BBAA <- dstats_BBAA %>%
    filter(p_adjusted <= 0.05) %>%   # remove trios without significant D
    mutate(pair = pmap_chr(., f_sorted_pair)) %>%
    group_by(pair) %>%
    mutate(d_is_max = f4_ratio == max(f4_ratio)) %>%
    filter(d_is_max) %>%
    mutate(f_is_max = f4_ratio == max(f4_ratio)) %>%
    filter(f_is_max) %>%
    mutate(p_is_min = p_adjusted == min(p_adjusted)) %>%
    filter(p_is_min) %>%
    ungroup()
)

(minmax_Dmin <- dstats_Dmin %>%
    filter(p_adjusted <= 0.05) %>%   # remove trios without significant D
    mutate(pair = pmap_chr(., f_sorted_pair)) %>%
    group_by(pair) %>%
    mutate(d_is_max = f4_ratio == max(f4_ratio)) %>%
    filter(d_is_max) %>%
    mutate(f_is_max = f4_ratio == max(f4_ratio)) %>%
    filter(f_is_max) %>%
    mutate(p_is_min = p_adjusted == min(p_adjusted)) %>%
    filter(p_is_min) %>%
    ungroup()
)


### Topology overlap and proportions
(sign_BBAA <- dstats_BBAA %>% 
  filter(p_adjusted < 0.05))

(sign_Dmin <- dstats_Dmin %>% 
  filter(p_adjusted < 0.05))

semi_join(sign_BBAA[1:3], sign_Dmin[1:3])   # intersection (significant trios in common)
full_join(sign_BBAA[1:3], sign_Dmin[1:3])   # all significant trios (52 %)

trios <- combinations(n = length(populations), r = 3, v = populations)
duos <- combinations(n = length(populations), r = 2, v = populations)

nrow(sign_BBAA) / nrow(trios)
nrow(sign_Dmin) / nrow(trios)
# nrow(bbaa_pairs) / nrow(duos)
# nrow(dmin_pairs) / nrow(duos)


### Convert to symmetrical
(fsym_BBAA <- minmax_BBAA %>%
   select(P2, P3, f4_ratio, p_adjusted) %>%
   bind_rows(minmax_BBAA %>% select(P2 = P3, P3 = P2, f4_ratio, p_adjusted)) %>%
   arrange(desc(f4_ratio))
)

(fsym_Dmin <- minmax_Dmin %>%
    select(P2, P3, f4_ratio, p_adjusted) %>%
    bind_rows(minmax_Dmin %>% select(P2 = P3, P3 = P2, f4_ratio, p_adjusted)) %>%
    arrange(desc(f4_ratio))
)


### Create tile plot (f4 ratio)
(fb <- fsym_BBAA %>%
    ggplot(aes(x = fct_relevel(P2, pop_by_region), 
               y = fct_relevel(P3, pop_by_region))) +
    geom_tile(aes(fill = f4_ratio), color = "white") +
    labs(title = NULL,
         x = NULL,
         y = NULL) +
    scale_fill_viridis_c(option = "plasma",
                         direction = -1,
                         limits = c(min(fsym_BBAA$f4_ratio), max(fsym_BBAA$f4_ratio)),
                         na.value = "white",
                         name = "f4 ratio\nBBAA topology") +
    scale_x_discrete(position = "top", labels = function(x) str_replace(str_replace(x, "tan", "sp1"), "esp", "sp2")) +
    scale_y_discrete(limits = rev, labels = function(y) str_replace(str_replace(y, "tan", "sp1"), "esp", "sp2")) +
    coord_equal(clip = "off") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(size = 2.75, hjust = 0, angle = 45, family="mono"),
          axis.text.y = element_text(size = 2.75, hjust = 0, family="mono"),
          legend.title = element_text(face = "bold", vjust = 5),
          legend.key.size = unit(4, "mm"),
          panel.border = element_rect(colour = "grey20", fill = NA, linewidth = 0.25),
          panel.grid = element_line(colour = "grey90", linewidth = 0.1),
          plot.margin =  margin(l = 2, b = 15),
          legend.margin = margin(l = 5, r = 2)
    )
)


(fd <- fsym_Dmin %>%
    ggplot(aes(x = fct_relevel(P2, pop_by_region), 
               y = fct_relevel(P3, pop_by_region))) +
    geom_tile(aes(fill = f4_ratio), color = "white") +
    labs(title = NULL,
         x = NULL,
         y = NULL) +
    scale_fill_viridis_c(option = "plasma",
                         direction = -1,
                         limits = c(min(fsym_BBAA$f4_ratio), max(fsym_BBAA$f4_ratio)),
                         na.value = "white",
                         name = "f4 ratio\nDmin topology") +
    scale_x_discrete(position = "top", 
                     labels = function(x) str_replace(str_replace(x, "tan", "sp1"), "esp", "sp2")) +
    scale_y_discrete(limits = rev, 
                     labels = function(y) str_replace(str_replace(y, "tan", "sp1"), "esp", "sp2")) +
    coord_equal(clip = "off") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(size = 2.75, hjust = 0, angle = 45, family="mono"),
          axis.text.y = element_text(size = 2.75, hjust = 0, family="mono"),
          legend.title = element_text(face = "bold", vjust = 5),
          legend.key.size = unit(4, "mm"),
          panel.border = element_rect(colour = "grey20", fill = NA, linewidth = 0.25),
          panel.grid = element_line(colour = "grey90", linewidth = 0.1),
          plot.margin =  margin(l = 2, b = 15),
          legend.margin = margin(l = 5, r = 2)
    )
)


### Combine panels
(f_comb <- fb /
    fd  +
    plot_annotation(tag_levels = c("A")) &
    theme(plot.tag = element_text(face = "bold"))
)


### Save heatmap
ggsave(
  filename = "figures/FigS7_f4ratio.pdf",
  plot = f_comb,
  width = 10,
  height = 15,
  units = "cm",
  device = cairo_pdf
)



### ===========================================================================================
### Alternative figure versions (Patterson's D)

### Convert to symmetrical
# (dsym_BBAA <- minmax_BBAA %>%
#    select(P2, P3, Dstatistic, p_adjusted) %>%
#    bind_rows(minmax_BBAA %>% select(P2 = P3, P3 = P2, Dstatistic, p_adjusted)) %>%
#    arrange(desc(Dstatistic))
# )
# 
# (dsym_Dmin <- minmax_Dmin %>%
#     select(P2, P3, Dstatistic, p_adjusted) %>%
#     bind_rows(minmax_Dmin %>% select(P2 = P3, P3 = P2, Dstatistic, p_adjusted)) %>%
#     arrange(desc(Dstatistic))
# )
# 
# 
# ### Create tile plot (D)
# (db <- dsym_BBAA %>%
#     ggplot(aes(x = fct_relevel(P2, pop_by_region), 
#                y = fct_relevel(P3, pop_by_region))) +
#     geom_tile(aes(fill = Dstatistic), color = "white") +
#     labs(title = NULL,
#          x = NULL,
#          y = NULL) +
#     scale_fill_viridis_c(option = "plasma",
#                          direction = -1,
#                          limits = c(min(dsym_BBAA$Dstatistic), max(dsym_BBAA$Dstatistic)),
#                          na.value = "white",
#                          name = "D-statistic\nBBAA topology") +
#     scale_x_discrete(position = "top", 
#                      labels = function(x) str_replace(str_replace(x, "tan", "sp1"), "esp", "sp2")) +
#     scale_y_discrete(limits = rev, 
#                      labels = function(y) str_replace(str_replace(y, "tan", "sp1"), "esp", "sp2")) +
#     coord_equal(clip = "off") +
#     theme_minimal(base_size = 8) +
#     theme(axis.text.x = element_text(size = 3, hjust = 0, angle = 45),
#           axis.text.y = element_text(size = 3, hjust = 0),
#           legend.title = element_text(face = "bold", vjust = 5),
#           legend.key.size = unit(4, "mm"),
#           panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#           plot.margin =  margin(l = 2, b = 15),
#           legend.margin = margin(l = 5, r = 2)
#     )
# )
# 
# (dd <- dsym_Dmin %>%
#     ggplot(aes(x = fct_relevel(P2, pop_by_region), 
#                y = fct_relevel(P3, pop_by_region))) +
#     geom_tile(aes(fill = Dstatistic), color = "white") +
#     labs(title = NULL,
#          x = NULL,
#          y = NULL) +
#     scale_fill_viridis_c(option = "plasma",
#                          direction = -1,
#                          limits = c(min(dsym_BBAA$Dstatistic), max(dsym_BBAA$Dstatistic)),
#                          na.value = "white",
#                          name = "D-statistic\nDmin topology") +
#     scale_x_discrete(position = "top", 
#                      labels = function(x) str_replace(str_replace(x, "tan", "sp1"), "esp", "sp2")) +
#     scale_y_discrete(limits = rev, 
#                      labels = function(y) str_replace(str_replace(y, "tan", "sp1"), "esp", "sp2")) +
#     coord_equal(clip = "off") +
#     theme_minimal(base_size = 8) +
#     theme(axis.text.x = element_text(size = 3, hjust = 0, angle = 45),
#           axis.text.y = element_text(size = 3, hjust = 0),
#           legend.title = element_text(face = "bold", vjust = 5),
#           legend.key.size = unit(4, "mm"),
#           panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#           plot.margin =  margin(l = 2, b = 15),
#           legend.margin = margin(l = 5, r = 2)
#     )
# )
# 
# 
# ### Combine panels
# (d_comb <- db /
#     dd  +
#     plot_annotation(tag_levels = c("a"))
# )
# 
# 
# ### Save heatmap
# ggsave(
#   filename = "figures/FigS7_D.pdf",
#   plot = f_comb,
#   width = 10,
#   height = 15,
#   units = "cm",
#   device = cairo_pdf
# )
