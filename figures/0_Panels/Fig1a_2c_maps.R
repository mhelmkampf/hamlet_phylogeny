### ============================================================================
### Hamlet phylogeny: R code to plot Fig. 1a, sample map
### Created 21 Nov 2023 by Martin Helmkampf
### Modified by Floriane Coulmance for final figures
### usage:
### Rscript Fig1a_2c_maps.R
### ============================================================================


### Clear workspace
dev.off()
rm(list = ls())


### Load needed library
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(data.table)
library(scatterpie)
library(ggtext)
library(tidyverse)


### Load World map
world <- ne_countries(scale = "large", returnclass = "sf")


### Read in sample list
samples <- read_lines("/hamlet_phylogeny/metadata/samples_ids.txt") %>%
           as_tibble() %>% rename(Sample = value) %>%
           mutate(SampleID = str_sub(Sample, end = -7),
                  spec = str_sub(Sample, -6, -4),
                  loc = str_sub(Sample, -3, -1))


### Add coordinates to sample table
samples_plus <- samples %>%
                filter(!(spec %in% c("tor", "tab", "tig"))) %>%
                mutate(coord_N = case_when(loc == "ala" ~ 22.399490,
                                           loc == "arc" ~ 20.203940,
                                           loc == "are" ~ 22.115278,
                                           loc == "bar" ~ 13.223694,
                                           loc == "bel" ~ 16.765278,
                                           loc == "boc" ~ 9.332778,
                                           loc == "flk" ~ 24.752580,
                                           loc == "gun" ~ 9.290722,
                                           loc == "hai" ~ 19.677000,
                                           loc == "hon" ~ 16.030000,
                                           loc == "liz" ~ 19.155547,
                                           loc == "pri" ~ 17.952820,
                                           loc == "qui" ~ 20.978806,
                                           loc == "san" ~ 12.501617,
                                           loc == "tam" ~ 21.475878),
                       coord_W = case_when(loc == "ala" ~ -89.656060,
                                           loc == "arc" ~ -91.970770,
                                           loc == "are" ~ -91.398333,
                                           loc == "bar" ~ -59.64506,
                                           loc == "bel" ~ -88.14417,
                                           loc == "boc" ~ -82.25472,
                                           loc == "flk" ~ -80.76065,
                                           loc == "gun" ~ -78.14039,
                                           loc == "hai" ~ -71.843000,
                                           loc == "hon" ~ -83.32861,
                                           loc == "liz" ~ -95.863903,
                                           loc == "pri" ~ -67.05643,
                                           loc == "qui" ~ -86.800111,
                                           loc == "san" ~ -81.731367,
                                           loc == "tam" ~ -97.227000))


### Set species colors
colstandard <- read_tsv("/hamlet_phylogeny/metadata/species_colors.tsv", 
                        col_types = "cc", 
                        col_names = TRUE)
colstandard <- colstandard[order(colstandard$Species),]


### Create pie chart data
pie_dat <- dcast(samples_plus, loc ~ spec, fill = 0, value.var = "spec", fun.aggregate = length) %>%
           mutate(coord_N = case_when(loc == "ala" ~ 22.399490,
                                      loc == "arc" ~ 20.203940,
                                      loc == "are" ~ 22.115278,
                                      loc == "bar" ~ 13.223694,
                                      loc == "bel" ~ 16.765278,
                                      loc == "boc" ~ 9.332778,
                                      loc == "flk" ~ 24.752580,
                                      loc == "gun" ~ 9.290722,
                                      loc == "hai" ~ 19.677000,
                                      loc == "hon" ~ 16.030000,
                                      loc == "liz" ~ 19.155547,
                                      loc == "pri" ~ 17.952820,
                                      loc == "qui" ~ 20.978806,
                                      loc == "san" ~ 12.501617,
                                      loc == "tam" ~ 21.475878),
                  coord_W = case_when(loc == "ala" ~ -89.656060,
                                      loc == "arc" ~ -91.970770,
                                      loc == "are" ~ -91.398333,
                                      loc == "bar" ~ -59.64506,
                                      loc == "bel" ~ -88.14417,
                                      loc == "boc" ~ -82.25472,
                                      loc == "flk" ~ -80.76065,
                                      loc == "gun" ~ -78.14039,
                                      loc == "hai" ~ -71.843000,
                                      loc == "hon" ~ -83.32861,
                                      loc == "liz" ~ -95.863903,
                                      loc == "pri" ~ -67.05643,
                                      loc == "qui" ~ -86.800111,
                                      loc == "san" ~ -81.731367,
                                      loc == "tam" ~ -97.227000),
                  Atlantic = case_when(loc == "ala" ~ 0,
                                       loc == "arc" ~ 0,
                                       loc == "are" ~ 0,
                                       loc == "bar" ~ 11,
                                       loc == "bel" ~ 0,
                                       loc == "boc" ~ 0,
                                       loc == "flk" ~ 0,
                                       loc == "gun" ~ 0,
                                       loc == "hai" ~ 3,
                                       loc == "hon" ~ 0,
                                       loc == "liz" ~ 0,
                                       loc == "pri" ~ 40,
                                       loc == "qui" ~ 0,
                                       loc == "san" ~ 0,
                                       loc == "tam" ~ 0),
                  GOM = case_when(loc == "ala" ~ 2,
                                  loc == "arc" ~ 20,
                                  loc == "are" ~ 3,
                                  loc == "bar" ~ 0,
                                  loc == "bel" ~ 0,
                                  loc == "boc" ~ 0,
                                  loc == "flk" ~ 19,
                                  loc == "gun" ~ 0,
                                  loc == "hai" ~ 0,
                                  loc == "hon" ~ 0,
                                  loc == "liz" ~ 6,
                                  loc == "pri" ~ 0,
                                  loc == "qui" ~ 0,
                                  loc == "san" ~ 0,
                                  loc == "tam" ~ 4),
                  Caribbean = case_when(loc == "ala" ~ 0,
                                        loc == "arc" ~ 0,
                                        loc == "are" ~ 0,
                                        loc == "bar" ~ 0,
                                        loc == "bel" ~ 66,
                                        loc == "boc" ~ 46,
                                        loc == "flk" ~ 0,
                                        loc == "gun" ~ 16,
                                        loc == "hai" ~ 0,
                                        loc == "hon" ~ 75,
                                        loc == "liz" ~ 0,
                                        loc == "pri" ~ 0,
                                        loc == "qui" ~ 12,
                                        loc == "san" ~ 4,
                                        loc == "tam" ~ 0),
                  sum = as.numeric(rowSums(.[, 2:21])),
                  radius = log(sum) / 1.8)


### Plot map
(c <- world %>% ggplot() +
    geom_sf(color = "grey20", fill = "grey90") +
      ## define limits of the map to span central america
    coord_sf(xlim = c(-100, -50), 
             ylim = c(5, 30), 
             expand = FALSE) +
      ## axis labels
    xlab("Longitude") + ylab("Latitude") +
      ## pie charts with parts as species count
      ## uncomment the next lines for Fig1a
    # geom_scatterpie(aes(x = coord_W, y = coord_N, group = loc, r = radius),
    #                 data = pie_dat, cols = colnames(pie_dat)[2:21], color = NA, alpha = 1) +
      ## pie charts as biogeographic regions
      ## uncomment the next lines for Fig2ca
    geom_scatterpie(aes(x = coord_W, y = coord_N, group = loc, r = radius),
                    data = pie_dat, cols = colnames(pie_dat)[24:26], color = NA, alpha = 1) +
      ## pie chart legend
      ## uncomment the next lines for Fig1a
    # geom_scatterpie_legend(pie_dat$radius, x = -56, y = 26, n = 3, color = "black", breaks = c(min(pie_dat$radius), mean(pie_dat$radius), max(pie_dat$radius)),
    # labeller = function(x) round(exp(x * 1.8), digits = 0)) +
      ## color pie charts by species
      ## uncomment the next lines for Fig1a
    # scale_fill_manual(values = colstandard$Color,
    #                   labels = c("abe" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_aberrans.l.cairo.png' width='55' /><br>*H. aberrans*",
    #                              "aff" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_affinis.l.cairo.png' width='55' /><br>*H. affinis*",
    #                              "atl" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_atlahua.l.cairo.png' width='55' /><br>*H. atlahua*",
    #                              "cas" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_castroaguirrei.l.cairo.png' width='55' /><br>*H. castroaguirrei*",
    #                              "chl" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_chlorurus.l.cairo.png' width='55' /><br>*H. chlorurus*",
    #                              "eco" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_ecosur.l.cairo.png' width='55' /><br>*H. ecosur*",
    #                              "esp" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_espinoza.l.cairo.png' width='55' /><br>*H. sp2*",
    #                              "flo" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_floridae.l.cairo.png' width='55' /><br>*H. floridae*",
    #                              "gem" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_gemma.l.cairo.png' width='55' /><br>*H. gemma*",
    #                              "gum" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_gumigutta.l.cairo.png' width='55' /><br>*H. gummigutta*",
    #                              "gut" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_guttavarius.l.cairo.png' width='55' /><br>*H. guttavarius*",
    #                              "ind" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_indigo.l.cairo.png' width='55' /><br>*H. indigo*",
    #                              "lib" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_liberte.l.cairo.png' width='55' /><br>*H. liberte*",
    #                              "may" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_maya.l.cairo.png' width='55' /><br>*H. maya*",
    #                              "nig" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_nigricans.l.cairo.png' width='55' /><br>*H. nigricans*",
    #                              "pro" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_providencianus.l.cairo.png' width='55' /><br>*H. providencianus*",
    #                              "pue" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_puella.l.cairo.png' width='55' /><br>*H. puella*",
    #                              "ran" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_randallorum.l.cairo.png' width='55' /><br>*H. randallorum*",
    #                              "tan" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_sp.l.cairo.png' width='55' /><br>*H. sp1*",
    #                              "uni" = "<img src='/hamlet_phylogeny/metadata/hamlet_logos/H_unicolor.l.cairo.png' width='55' /><br>*H. unicolor*")) +
    ## color pie charts by biogeographic regions
    ## uncomment for Fig2c
  scale_fill_manual(values = c("olivedrab4", "coral2", "royalblue1"),
                    labels = c("Atlantic", "Gulf of Mexico", "Caribbean")) +
    ## add names of different seas
  annotate(geom = "text",
             x = c(-90, -74, -62),
             y = c(26, 14.5, 23),
             label = c("Gulf of Mexico", "Caribbean", "Atlantic"),
             fontface = "italic",
             color = "grey50", size = 4) +
      ## add name acronyms of all sites near piecharts
    annotate(geom = "text",
             x = c(-88.656,
                   -93.871,
                   -91.398,
                   -57.145,
                   -90.844,
                   -82.255,
                   -83.261,
                   -76.140,
                   -70.843,
                   -79.629,
                   -95.864,
                   -65.056,
                   -84.300,
                   -79.731,
                   -96.227),
             y = c(23.399,
                   22.104,
                   23.415,
                   13.224,
                   14.765,
                   6.6328,
                   24.753,
                   11.291,
                   20.677,
                   16.030,
                   17.456,
                   15.953,
                   20.979,
                   12.502,
                   22.876),
             label = c("ala",
                       "arc",
                       "are",
                       "bar",
                       "bel",
                       "boc",
                       "flk",
                       "gun",
                       "hai",
                       "hon",
                       "liz",
                       "pri",
                       "qui",
                       "san",
                       "tam"),
             fontface = "italic",
             color = "grey20", size = 3) +
      ## arrange species name legend in 2 rows
    guides(fill = guide_legend(nrow = 2)) +
      ## adjust overall esthetics
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 10),
          legend.text =  element_markdown(size = 10),
          legend.key = element_blank(),
          legend.box = "vertical",
          legend.position = "bottom",
          legend.title = element_blank())
)


### Save as PDF
### uncomment for Fig1a
# ggsave(
#   filename = "/hamlet_phylogeny/figures/Fig1a_map.pdf",
#   plot = c,
#   width = 36,
#   height = 21.5,
#   units = "cm",
#   device = cairo_pdf
# )


### Save as PDF
### uncomment for Fig2c
ggsave(
  filename = "/hamlet_phylogeny/figures/Fig2c_legend.pdf",
  plot = c,
  width = 36,
  height = 21.5,
  units = "cm",
  device = cairo_pdf
)
