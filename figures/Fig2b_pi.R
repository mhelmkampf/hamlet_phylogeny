### ============================================================================
### Hamlet phylogeny: R code to plot Fig. 2b (nucleotide diversity, tskit)
### Created Jan 2024 by Martin Helmkampf
### ============================================================================


### Preparations
library(tidyverse)
library(matrixStats)
library(ggpattern)
library(cowplot)

# Set path to root directory of git repository "hamlet_phylogeny"

rm(list = ls())


### Define dataset
ds <- "phylo2e"
level <- "pop"
lgs <- paste0("LG", sprintf('%0.2d', 1:24))


### Declare populations
samples <- read_tsv("metadata/samples_ids.txt", col_names = "Sample") %>%
  mutate(Population = str_sub(Sample, -6, -1))
# samples <- samples %>% filter(!grepl("tor|tab|tig", Population))   # execute if outgroups were removed

populations <- levels(factor(samples$Population))


### Define regions
gulf <- c("liz", "tam", "ala", "arc", "are", "flk")
caribbean <- c("bel", "boc", "gun", "hon", "san", "qui")
atlantic <- c("bar", "hai", "pri")


### Read in data from files
df <- data.frame(0)

for (region in lgs) {
  pi <- as.numeric(read_lines(file = paste0("results/pi/pi_", ds, "_", region, "_", level, "_n1.txt")))
  df <- cbind(df, pi)
}

df <- subset(df, select = -1)

colnames(df) <- lgs


### Calculate means and sd, add population and region ids
nucdiv <- as_tibble(df) %>%
  mutate(mean_pi = rowMeans(.),
         sd_pi = rowSds(as.matrix(.))) %>%
  select(mean_pi, sd_pi) %>%
  cbind(populations) %>%
  rename(population = populations) %>%
  mutate(spec = str_sub(population, -6, -4),
         loc = str_sub(population, -3, -1),
         region = case_when(
           loc %in% gulf ~ "Gulf of Mexico",
           loc %in% caribbean ~ "Caribbean",
           loc %in% atlantic ~ "Atlantic",
           TRUE ~ "ungrouped"),
         clade = case_when(
           spec %in% c("atl", "eco", "flo") ~ "Small clade",
           TRUE ~ "Large clade"
         )
  ) %>%
  filter(!grepl("tor|tab|tig", population),
         !population %in% samples)


### Create histogram plot
(g <- nucdiv %>%
    ggplot(aes(x = reorder(population, -mean_pi),
               y = mean_pi,
               fill = region,
               pattern = clade)) +
    geom_bar_pattern(position = "dodge",
             stat = "identity",
             color = "grey20",
             alpha = 0.8,
             lwd = 0.4,
             width = 1,
             pattern_color = "grey35",
             pattern_fill = "grey35",
             pattern_angle = 45,
             pattern_density = 0.02,
             pattern_spacing = 0.02,
             pattern_key_scale_factor = 0.75) +
    geom_errorbar(aes(ymin = mean_pi - sd_pi, ymax = mean_pi + sd_pi),
                  color = "grey50",
                  width = 0.4) +
    scale_fill_manual(breaks = c("Gulf of Mexico", "Caribbean", "Atlantic"),
                      values = c("coral2", "royalblue1", "olivedrab4")) +
    scale_pattern_manual(breaks = c("Small clade", "Large clade"), 
                         values = c("stripe", "none")) +
    labs(# title = paste0("Mean genome-wide nucleotide diversity (π) per population"), 
         x = NULL, 
         y = expression(paste("Nucleotide diversity (", italic("\u03c0"), ")"))
         ) +
    guides(fill = "none",
           pattern = guide_legend(title = "Clade", override.aes = list(fill = "white"))
           ) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12, margin = margin(0, 0, 5, 0)),
          # axis.text.x = element_text(hjust = 1, vjust = 1, angle = 35, size = 7,
          #                            margin = unit(c(-3, 0, 0, 0), "mm")),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16, vjust = 3),
          legend.position = c(0.78, 0.875), legend.box = "horizontal", legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          plot.margin = unit(c(3, 3, 0, 3), "mm")
          )
)


### Set species colors
colstandard <- read_tsv("metadata/species_colors.tsv", 
                        col_types = "cc", 
                        col_names = TRUE)


### Create species color table
scol <- nucdiv %>%
  arrange(-mean_pi) %>%
  mutate(position = seq(1:length(nucdiv$population)))


### Plot species colors
(s <- ggplot(scol, aes(xmin = position, xmax = position + 1, 
                       ymin = -Inf, ymax = Inf)) +
  geom_rect(aes(fill = spec)) +
    scale_x_discrete(breaks = seq(1:length(nucdiv$population))) +
    scale_y_continuous(breaks = .5, limits = c(0, 1)) +
    scale_fill_manual(values = colstandard$Color) +
    labs(x = NULL, 
         y = "Species") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          # axis.title.y = element_text(angle = 0, vjust = 0.5, 
          #                             color = "gray20", size = 9,
          #                             margin = margin(r = -40)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(-3, 0, 3, 0), "mm")
    )
)


### Combine plots
(p_comb <- plot_grid(g, s, 
                     align = "v", 
                     nrow = 2, 
                     rel_heights = c(0.93, 0.07),
                     scale = c(1, 1.015))
)


### Save plot
ggsave(
  filename = paste0("figures/Fig2b_pi.png"),
  plot = p_comb,
  width = 20,
  height = 10,
  units = "cm",
  device = png,
  bg = "white"
)
