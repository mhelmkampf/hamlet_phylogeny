### ============================================================================
### Hamlet phylogeny: R code to plot Fig. S9 (Admixture analysis k=2 to k=8)
### Created Nov 2023 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)
library(cowplot)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Set colors
colstandard <- read_tsv("metadata/species_colors.tsv",
                        col_types = "cc",
                        col_names = TRUE)


### Define species and order
gulf_sp <- c("atl", "cas", "eco", "esp", "flo", "gem")
caribbean_sp <- c("aff", "chl", "gum", "gut", "ind", "lib", "may", "pro", "tan")
both_sp <- c("abe", "nig", "pue", "uni", "ran")

order_sp <- c(gulf_sp, both_sp, caribbean_sp)


### Read in and process data for k = 2:8
admix_data <- map(2:8, function(k) {
  read_delim(paste0("results/admix/AdmcvProp_phyps2e_m2n1l5_k", k, ".tsv"), delim = " ", col_names = "Sample") %>%
    mutate(Species = str_sub(Sample, -6, -4),
           Location = str_sub(Sample, -3, -1)) %>%
    pivot_longer(cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion") %>%
    mutate(k = k)
})


### Generate ordered sample list
samples <- admix_data[[1]] %>% distinct(Sample, Species, Location)
samples$Species <- factor(samples$Species, levels = order_sp)
ord_sp <- samples[order(samples$Species, samples$Location), ] %>% pull(Sample)


### Plotting function
plot_admixture <- function(data, k, colors = NULL) {
  ancestries <- unique(data$Ancestry)
  
  if (is.null(colors)) {
    ncolors <- length(ancestries)
    colors <- colorRampPalette(brewer.pal(9, "Set1"))(ncolors)
  }
  
  names(colors) <- ancestries
  
  ggplot(data, aes(x = fct_relevel(Sample, ord_sp), y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = NULL, tag = paste0("k = ", k)) +
    scale_fill_manual(values = colors) +
    theme_minimal(base_size = 11) +
    theme(text = element_text(color = "grey20"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5),
          legend.position = "none",
          plot.margin = unit(c(1, 10, 1, 1), "mm"))
}


### Manually define color palettes
admix_colors <- list(
  `2` = c("X1" = "#377EB8", "X2" = "#E41A1C"),
  `3` = c("X1" = "#4DAF4A", "X2" = "#377EB8", "X3" = "#E41A1C"),
  `4` = c("X1" = "#4DAF4A", "X2" = "#377EB8", "X3" = "#E41A1C", "X4" = "#984EA3"),
  `5` = c("X1" = "#984EA3", "X2" = "#FF7F00", "X3" = "#4DAF4A", "X4" = "#E41A1C", "X5" = "#377EB8"),
  `6` = c("X1" = "#984EA3", "X2" = "#E41A1C", "X3" = "#4DAF4A", "X4" = "#377EB8", "X5" = "#FF7F00", "X6" = "#FFFF33"),
  `7` = c("X1" = "#E41A1C", "X2" = "#4DAF4A", "X3" = "#984EA3", "X4" = "#377EB8", "X5" = "#FF7F00", "X6" = "#A65628", "X7" = "#FFFF33"),
  `8` = c("X1" = "#377EB8", "X2" = "#A65628", "X3" = "#999999", "X4" = "#4DAF4A", "X5" = "#E41A1C", "X6" = "#FFFF33", "X7" = "#FF7F00", "X8" = "#984EA3")
)

### Store plots in list
plots <- map2(admix_data, 2:8, ~plot_admixture(.x, k = .y, colors = admix_colors[[as.character(.y)]]))


### Create species bar plot
coord_sp <- samples[order(samples$Species, samples$Location), ] %>%
  mutate(Ord_nr = row_number()) %>%
  select(Sample, Species, Ord_nr)

blocks_sp <- coord_sp %>%
  group_by(Species) %>%
  summarise(Start = min(Ord_nr) - 1, End = max(Ord_nr)) %>%
  left_join(colstandard) %>%
  replace(is.na(.), "gray80")

(p_sp <- ggplot(blocks_sp, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf)) +
  geom_rect(aes(fill = Species), color = "grey20", lwd = 0) +
  scale_y_continuous(breaks = .5, limits = c(0, 1)) +
  scale_x_discrete(breaks = coord_sp$Ord_nr,
                   labels = coord_sp$Species, expand = c(0, 0)) +
  scale_fill_manual(values = colstandard$Color) +
  labs(tag = "Species") +
  theme_minimal(base_size = 11) +
  theme(text = element_text(color = "grey20"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        plot.tag = element_text(angle = -90),
        plot.tag.position = c(1.02, 0.5),
        plot.margin = unit(c(4, 10, 3, 1), "mm"))
)



### ============================================================================
### Plot cross-validation error

cv_error <- data.frame(
  k = 1:12,
  Error = c(0.59175, 
            0.58147, 
            0.57689, 
            0.57688, 
            0.57609, 
            0.57783, 
            0.58136, 
            0.58323, 
            0.58668, 
            0.59434, 
            0.59585, 
            0.61051)
)

# Plot with ggplot2
(p_cv <- ggplot(cv_error, aes(x = k, y = Error)) +
    geom_line(color = "grey40", size = 1) +
    geom_point(color = "grey40", size = 2) +
    scale_x_continuous(breaks = 1:12) +
    labs(x = "k",
         y = NULL,
         tag = "CV error") +
    theme_minimal(base_size = 11) +
    theme(text = element_text(color = "grey20"),
          panel.grid.minor = element_blank(),
          plot.tag = element_text(angle = -90),
          plot.tag.position = c(1.02, 0.5))
)


### Combine all plots
(p_comb <- plot_grid(
  plots[[1]], plots[[2]], plots[[3]], plots[[4]],  plots[[5]], plots[[6]], plots[[7]],
  p_sp,
  p_cv,
  align = "v",
  nrow = 9,
  rel_heights = c(rep(1, 7), 0.5, 3)
))


### Save as PDF
ggsave("figures/FigS9_admixk2-8.pdf", plot = p_comb, width = 20, height = 26, units = "cm", device = cairo_pdf)

### Save as PNG
ggsave("figures/FigS9_admixk2-8.png", plot = p_comb, width = 20, height = 26, units = "cm", device = png)



### ============================================================================
### Region bar (not functional)

### Define regions
# gulf <- c("liz", "tam", "ala", "arc", "are", "flk")
# west_carib <- c("bel", "boc", "gun", "hon", "san", "qui")
# east_carib <- c("bar", "hai", "pri")
# 
# region_colors <- tibble(
#   Region = c("Eastern Caribbean", "Western Caribbean", "Gulf of Mexico"),
#   Color = c("olivedrab4", "royalblue1", "coral2")
# )

### Create region bar plot
# coord_rg <- samples[order(samples$Species, samples$Region), ] %>%
#   mutate(Ord_nr = row_number()) %>%
#   select(Sample, Region, Ord_nr)
# 
# blocks_rg <- coord_rg %>%
#   left_join(region_colors) %>%
#   replace(is.na(.), "gray80")
# 
# (p_sp <- ggplot(blocks_rg, aes(xmin = Ord_nr, xmax = Ord_nr+1, ymin = -Inf, ymax = Inf)) +
#     geom_rect(aes(fill = Region), color = "grey20", lwd = 0) +
#     scale_y_continuous(breaks = .5, limits = c(0, 1)) +
#     scale_x_discrete(breaks = coord_rg$Ord_nr,
#                      labels = coord_rg$Region, expand = c(0, 0)) +
#     scale_fill_manual(values = blocks_rg$Color) +
#     theme_minimal() +
#     theme(panel.grid = element_blank(),
#           axis.title = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           legend.position = "none",
#           plot.margin = unit(c(1, 10, 3, 1), "mm"))
# )
