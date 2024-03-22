### ============================================================================
### Hamlet phylogeny: R code to plot Fig. S1 (sample overview)
### Created Oct 2023 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)
library(grid)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Read in sample list
samples <- read_tsv(file = "metadata/samples_phylo.ids", col_names = "Sample")


### Relabel samples (samples file has been updated, keep for documentation purposes)
# samples$Sample <- samples$Sample %>%
#   str_replace(pattern = "54761nigliz", "54761atlliz") %>%
#   str_replace(pattern = "52989nigtam", "52989atltam") %>%
#   str_replace(pattern = "27936flotam", "27936atltam") %>%
#   str_replace(pattern = "PL17_35puepri", "PL17_35indpri") %>%
#   str_replace(pattern = "19294tangun", "19294affgun") %>%
#   str_replace(pattern = "62585uniarc", "62585esparc") %>%
#   str_replace(pattern = "62571uniarc", "62571esparc") %>%
#   str_replace(pattern = "62555uniarc", "62555esparc") %>%
#   str_replace(pattern = "27698gutqui", "27698abequi")
# 
# metaseq$Sample <- metaseq$Sample %>%
#   str_replace(pattern = "54761nigliz", "54761atlliz") %>%
#   str_replace(pattern = "52989nigtam", "52989atltam") %>%
#   str_replace(pattern = "27936flotam", "27936atltam") %>%
#   str_replace(pattern = "PL17_35puepri", "PL17_35indpri") %>%
#   str_replace(pattern = "19294tangun", "19294affgun") %>%
#   str_replace(pattern = "62585uniarc", "62585esparc") %>%
#   str_replace(pattern = "62571uniarc", "62571esparc") %>%
#   str_replace(pattern = "62555uniarc", "62555esparc") %>%
#   str_replace(pattern = "27698gutqui", "27698abequi")


### Add species / location information
samplesSL <- samples %>%
  mutate(Species = str_sub(Sample, -6, -4),
         Location = str_sub(Sample, -3, -1),
         Population = paste0(Species, Location)
  )


### Define groups
ingroup <- c("abe", "aff", "atl", "cas", "chl", "eco", "esp", "flo", "gem", "gum", 
             "gut", "ind", "lib", "may", "nig", "pro", "pue", "ran", "tan", "uni")
outgroup <- c("tab", "tig", "tor")

gulf <- c("liz", "tam", "ala", "arc", "are", "flk")
caribbean <- c("bel", "boc", "gun", "hon", "san", "qui")
atlantic <- c("bar", "hai", "pri")


### Summarize
counts <- dplyr::count(samplesSL, Species, Location) %>%
  replace(is.na(.), 0) %>%
  arrange(Species)


### Order within groups
order_sp <- samplesSL %>%
  mutate(group_sp = case_when(
    Species %in% ingroup ~ "ingroup", 
    Species %in% outgroup ~ "outgroup")) %>%
  dplyr::count(Species, group_sp) %>%
  group_by(group_sp) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  pull(Species)

order_loc <- samplesSL %>%
  mutate(Region = case_when(
    Location %in% gulf ~ "Gulf of Mexico",
    Location %in% caribbean ~ "Caribbean",
    Location %in% atlantic ~ "Atlantic")) %>%
  dplyr::count(Location, Region) %>%
  group_by(Region) %>%
  arrange(match(Region, c("Caribbean", "Atlantic", "Gulf of Mexico")), desc(n)) %>%
  pull(Location)


### Set species colors and text objects
colstandard <- read_tsv("metadata/species_colors.tsv", 
                        col_types = "cc", 
                        col_names = TRUE)

scol <- left_join(tibble(Species = order_sp), colstandard) %>%
  replace(is.na(.), "gray80")

sum <- textGrob(sum(counts$n), gp = gpar(fontsize = 12, fontface = "bold"))


### Plot overview
(plot <- counts %>%
  mutate(
    Location = fct_relevel(Location, order_loc),
    Species = fct_relevel(Species, rev(order_sp))
  ) %>%
  ggplot(aes(x = Location, y = Species)) +
  geom_count(aes(size = n, color = Species), show.legend = FALSE) +
  geom_text(aes(label = n), size = 2.75, nudge_x = 0.46, color = "gray50") + # add counts as text
  geom_vline(xintercept = c(6.6, 9.6), col = "gray80", linetype = "dashed") +
  geom_hline(yintercept = 3.5, col = "gray80", linetype = "dashed") +
  scale_color_manual(values = rev(scol$Color)) +
  scale_size_area(max_size = 11) +
  coord_cartesian(clip = "off") +
  labs(title = NULL, x = NULL, y = NULL) +
  annotate(geom = "text", 
           x = c(3.5, 8, 12.5, 16, 16), 
           y = c(24, 24, 24, 12.5, 2), 
           label = c("Caribbean", "Atlantic", "Gulf of Mexico", "Ingroup", "Outgroup"),
           angle = c(0, 0, 0, 270, 270),
           color = "gray20") +
  annotation_custom(sum, xmin = 16, xmax = 16, ymin = 0.07, ymax = 0.07) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic"),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.25, l = 0.25, unit = "cm")
  )
)


### Save to file
ggsave(
  filename = "figures/FigS1_tally.pdf",
  plot = plot,
  width = 22,
  height = 22,
  units = "cm",
  device = cairo_pdf
)


### Create table
# tally_tmp <- counts %>%
#   spread(Location, n) %>%
#   replace(is.na(.), 0) %>%
#   mutate(sumSp = rowSums(.[2:ncol(.)])) %>%
#   arrange(factor(Species, levels = order_sp)) %>%
#   add_row(summarize_all(., ~ if (is.numeric(.)) { 
#                                  sum(.) }
#                              else "sumLoc"))  # add column sums (= by Location)
# 
# tally <- tally_tmp[, match(order_loc, colnames(tally_tmp))] %>%
#   add_column(Species = tally_tmp$Species, .before = .[[1]]) %>%
#   add_column(sumSp = tally_tmp$sumSp)
# 
# write_tsv(tally, "data/phylo2e_tally.tsv")
