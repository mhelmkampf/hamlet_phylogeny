#!/usr/bin/env Rscript
# by: Floriane Coulmance: 16/11/2023
# usage:
# Rscript Fig3a-b_S12-19_gwas.R <data_path> <figure_path>
# -------------------------------------------------------------------------------------------------------------------
# <data_path> in : $BASE_DIR/outputs/gxp_clades/large/
# <figure_path> in : $BASE_DIR/figures/gxp_clades/large/
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(tidyverse)
library(hypogen)
library(hypoimg)
library(GenomicOriginsScripts)
library(plyr)
library(dplyr)
library(ggrepel)



# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:8]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
# data_path <- "/Users/fco/Desktop/PhD/2_CHAPTER2/chapter2/"
figure_path <- as.character(args[2]) # Path to the figure folder
# figure_path <- "/Users/fco/Desktop/PhD/2_CHAPTER2/chapter2/"



# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


prep_file <- function(f,path) {
  
  # Function to open and prepare file's columns 
  # relative to the type of GWAS file (50k window-averaged or snp-by-snp)
  
  # Identify window averaged or snp-by-snp files
  if (grepl(".10k.1k.txt.gz", f, fixed = TRUE)) {
    pattern <- "10k"
  } else if (grepl(".GxP.txt.gz", f, fixed = TRUE)) {
    pattern <- "GxP"
  } else if (grepl(".50k.5k.txt.gz", f, fixed = TRUE)) {
    pattern <- "50k"
  }
  
  # Prepare file and needed columns
  if (pattern == "50k") {
    run_files <- "50k_lmm"
    
    print(run_files)
    
    d <- read.table(paste0(path,f), header=TRUE) %>% left_join(hypo_chrom_start) %>% 
      mutate(RUN = run_files,
             LOG_P = AVG_p_wald,
             GPOS = MID_POS + GSTART)
    
  } else if (pattern == "10k") {
    run_files <- "10k_lmm"
    
    print(run_files)
    
    d <- read.table(paste0(path,f), header=TRUE) %>% left_join(hypo_chrom_start) %>% 
      mutate(RUN = run_files,
             LOG_P = AVG_p_wald,
             GPOS = MID_POS + GSTART) 
    
  } else if (pattern == "GxP") {
    run_files <- "snp_lmm"
    
    print(run_files)
    
    d <- read.table(paste0(path,f), header=TRUE) %>% left_join(hypo_chrom_start) %>%
      mutate(RUN = run_files,
             LOG_P = -log(p_wald, base = 10),
             MID_POS = POS,
             BIN_START = POS,
             BIN_END = POS,
             GPOS = POS + GSTART,
             RANGE = paste(CHROM, POS, sep="_"))
    
  } else {
    print("wrong file")
  }
  
  return(d)
  
}


threshold_table <- function(f) {
  
  # Create a table with regions of high association signal to plot
  
  assoc <- thresh  %>% select(CHROM, BIN_START, BIN_END, LOG_P, RUN) %>%
    setNames(., nm = c('chrom', 'start', 'end', 'log_p', 'run') ) %>%
    left_join(.,hypogen::hypo_chrom_start, by = c('chrom' = 'CHROM')) %>%
    mutate(gpos = (start+end)/2 + GSTART) %>%
    group_by(chrom) %>%
    mutate(check = gpos > lag(gpos,default = 0) + 50000,
           id = cumsum(check),
           gid = str_c(chrom,'_',id)) %>%
    group_by(gid) %>%
    dplyr::summarise(chrom = chrom[1],
                     start = min(start),
                     end = max(end),
                     gstart = min(start)+GSTART[1],
                     gend = max(end)+GSTART[1]) %>%
    mutate(gpos = (gstart+gend)/2)
  
  return(assoc)  
  
}


custom_annoplot_flo <- function (searchLG, xrange, genes_of_interest = cool_genes, genes_of_sec_interest = cool_genes,
                                 anno_rown = 3, width = 0.1, gene_color = 'darkgray', start, end) {
  
  # Extract annotation information of specific LG range and plot the specific region
  
  # Import annotation data for a specific LG (chromosome) and range
  df_list <- hypogen::hypo_annotation_get(searchLG = searchLG, xrange = xrange,
                                          genes_of_interest = genes_of_interest,
                                          genes_of_sec_interest = genes_of_sec_interest,
                                          anno_rown = anno_rown)
  
  print (df_list[[1]]$labelx)

  # Plot the annotation
  ggplot2::ggplot() +
    # add outlier area
    ggchicklet:::geom_rrect(data = tibble::tibble(start = start, end = end),
                       aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                       fill = "#E69F00", alpha = 0.15, r = unit(0.25, 'npc')) +
                       # fill=rgb(1,1,1,.3),color = rgb(1,1,1,.9)) +
    # add gene direction if known
    ggplot2::geom_segment(data = (df_list[[1]] %>% dplyr::filter(strand %in% c("+", "-"))),
                          aes(x = ps, xend = pe, y = yl, yend = yl, group = Parent),
                          lwd = 2, arrow = arrow(length = unit(20,"pt"), type = "closed"),
                          color = clr_genes) +
    # add gene extent if direction is unknowns
    ggplot2::geom_segment(data = (df_list[[1]] %>% dplyr::filter(!strand %in% c("+", "-"))),
                          aes(x = ps, xend = pe, y = yl, yend = yl, group = Parent),
                          lwd = 2, color = clr_genes) +
    # add exon 
    ggplot2::geom_rect(data = df_list[[2]],
                       aes(xmin = ps, xmax = pe, ymin = yl - width*2, ymax = yl + width*2, group = Parent),
                       fill = alpha(gene_color,.8)) +
    # add gene label
    ggplot2::geom_text(data = df_list[[1]],
                       aes(x = labelx, label = gsub("hpv1g000000", ".", label), y = yl - width*6), size = 20, color = clr_genes)
  
  
}


plot_panel_anno_flo <- function(outlier_id, label, lg, start, end,...)  {
  
  # Plot an annotation for outlier region of interest
  
  # Create plot title
  ttle <- stringr::str_sub(outlier_id,1,4) #%>% stringr::str_c(.,' (',project_inv_case(label),')')
  
  # Create plot of region annotation
  p <- custom_annoplot_flo(searchLG = lg,
                           xrange = c(start,end),
                           genes_of_interest = cool_genes,
                           anno_rown = 6, start = start, end = end) +
    # layout x ayis
    ggplot2::scale_x_continuous(name = ttle,
                                position = 'top',
                                expand = c(0,0),
                                limits = c(start, end),
                                labels = ax_scl,
                                breaks = NULL) +
    # layout y ayis
    ggplot2::scale_y_continuous(name = expression(bolditalic(Genes)), expand = c(0,.4)) +
    # use same boundaries for all panels
    ggplot2::coord_cartesian(xlim = c(start-window_buffer, end+window_buffer)) +
    # special panel layout for annotation panel
    ggplot2::theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   legend.position = 'none',
                   axis.title.x = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.x = element_text(size = 30),
                   plot.margin = ggplot2::margin(t = 1, r = 3, b = 3, 
                                                 l = 3))
  
  # Use correct greek symbols in labels if needed
  if(outlier_id == 'LG17_1') {
    p$layers[[5]]$data$label <- p$layers[[5]]$data$label %>%
      stringr::str_replace(., 'alpha', "\u03B1") %>%
      stringr::str_replace(.,  'beta', "\u03B2")
  }
  
  p
  
}


plot_panel_gxp <- function (lg, start, end, ...) {
  
  # Create gxp plot for 50k windowed GWAS results
  
  ggplot2::ggplot() +
    # colors the area under the curve & between 2 x-axis coordinates
    geom_ribbon(data = gxp_data2 %>% filter(CHROM == lg, MID_POS >= start, MID_POS <= end),
                aes(x = MID_POS, y = LOG_P, ymin = 0, ymax = LOG_P/0.20), fill=gxp_clr, alpha=.1) +
    # add snp-by-snp association points
    geom_point(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25),
               aes(x = MID_POS, y = LOG_P),
               size = 4, stroke = 0.2, color = "darkgray") +
    # add average over 10kb sliding windows association curve
    geom_line(data = gxp_data2 %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25),
              aes(x = MID_POS, y = LOG_P/0.20, color = RUN),
              size = 3, color = gxp_clr) +
    # add x-axis coordinates indicating the region of high association and used for phylogenetic tree
    scale_x_continuous(name = lg, position = "bottom", breaks = c(start, end)) +
    # y-axis esthetics
    scale_y_continuous(name = expression(bolditalic(-log[10](p))), expand = c(0, 0), limits = c(0, 30), sec.axis = sec_axis(~ . * 0.20),breaks = scales::pretty_breaks(n = 4)) +
    # determine legend looks
    guides(color = guide_legend(keyheight = unit(3,"pt"), keywidth = unit(20, "pt"),
                                override.aes = list(size = 2))) +
    # optional
    scale_color_manual(name = "GxP Trait", values = c("black", "#E69F00", gxp_clr)) +
    # use same boundaries for all panels
    ggplot2::coord_cartesian(xlim = c(start-window_buffer, end+window_buffer)) +
    # general esthetics    
    ggplot2::theme(panel.border = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  legend.position = "none",
                   axis.title.y = ggplot2::element_text(angle = 90, size = 40), 
                   axis.line.y = ggplot2::element_line(size = plot_lwd*1.5), 
                   axis.ticks.y = ggplot2::element_line(size = plot_lwd*1.5), 
                   axis.title.x = ggplot2::element_text(size = 50),
                   axis.text.x = ggplot2::element_text(size = 50), 
                  axis.line.x = element_line(),
                  axis.ticks.x = element_line(),
                   axis.text.y = ggplot2::element_text(size = 50),
                   axis.line.y.right = element_line(color = gxp_clr),
                   axis.text.y.right = element_text(color = gxp_clr),
                  axis.ticks.y.right = element_line(color = gxp_clr),
                   plot.margin = ggplot2::margin(t = 1, r = 3, b = 3, 
                                                 l = 3), ...)

}


plot_curt <- function (outlier_id, outlier_nr, lg, start, end, text = FALSE, label, ...) {
  
  # Create zoom plots into GWAS peak regions for one particular LG, with all necessary 
  
  # Annotation pannel
  p_g <- plot_panel_anno_flo(lg = lg, outlier_id = outlier_id, label = label,
                             start = start, end = end, genes = cool_genes)
  
  # Pannel for GWAS peak zoomed and 50kb averaged windows
  p_gxp <- plot_panel_gxp(lg = lg, start = start, end = end)
  
  # Assemble all panels
  p_curtain <- cowplot::plot_grid(p_g + no_title(), p_gxp + no_title(),
                                    ncol = 1, align = "v", rel_heights = c(1, 0.6), axis="tblr")
  
  # Save each individual zoom plots
  hypo_save(filename = paste0(figure_path, "/", outlier_id, ".pdf"),
            plot = p_curtain,
            width = 45,
            height = 12)
  
  print(p_curtain)
  
}


plot_regions <- function(region, outlier_list, label, path) {
  
  # Create plot for regions with highest associations within each LG
  
  # Create the plot
  p_single <- region %>% filter(outlier_id %in% outlier_list) %>%
    left_join(label) %>%
    dplyr::mutate(outlier_nr = row_number()) %>% 
    pmap(plot_curt, cool_genes = cool_genes) %>%
    cowplot::plot_grid(plotlist = ., nrow = row, rel_heights = heights,
                       labels = letters[1:length(outlier_list)] %>%
                         project_case())
  
  p_single
  
}



# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

# 1) Create whole genome Manhattan plot of GWAS association over 50kb sliding windows
# -------------------------------------------------------------------------------------------------------------------

gem50lmm <- "large.lmm.50k.5k.txt.gz"
gxp_data <- prep_file(gem50lmm,data_path)

(p <- ggplot() +
     geom_hypo_LG() +
     geom_point(data = gxp_data, aes(x = GPOS, y = LOG_P), size = .9) +
     scale_fill_hypo_LG_bg() +
     scale_x_hypo_LG(name = "Linkage Groups", position = "top") +
     scale_y_continuous(name = expression(italic('-log(p-value)'))) +
     theme_hypo() +
     theme(legend.position = 'none',
            axis.title.x = element_text(size =40),
            axis.title.y = element_text(size =40),
            axis.text.x.top = element_text(colour = 'darkgray', size = 40),
            axis.text.y = element_text(colour = 'darkgray', size = 50),
           axis.text.x = element_text(colour = 'darkgray', size = 50),
            plot.margin = unit(c(0,0,0,0), "cm")) 
)

hypo_save(filename = paste0(figure_path,"/gxp_large_gemma_lmm_plots.pdf"),
          plot = p,
          width = 45,
          height = 3.5)


# 2) Create zoomed in Manhattan plot of the 8 interesting regions identified by GWAS
# -------------------------------------------------------------------------------------------------------------------

  #Load the snp-by-snp file and format it 
gem10lmm <- "large.lmm.10k.1k.txt.gz"
gxp_data2 <- prep_file(gem10lmm,data_path)

  #Load the snp-by-snp file and format it 
gemSNPlmm <- "large.lmm.GxP.txt.gz"
gxp_snp <- prep_file(gemSNPlmm,data_path)

  # Set the threshold to retrieve 8 interesting regions
threshold <- 1.5

  # Select the 8 regions where the association signal is the highest
thresh <- gxp_data[gxp_data[, "LOG_P"] >= threshold,]
print(thresh)

  # Create table with the 8 regions of interest
region_table <- threshold_table(thresh) %>% 
  setNames(., nm = c("outlier_id", "lg", "start", "end", "gstart", "gend", "gpos", "heatmap"))

print(region_table)

  # Write 8 region of interest table to folder
write.table(region_table, file = paste0(data_path, "/gxp_large_gemma_lmm_regions.txt"), quote=FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  # List the 8 outlier regions' names
outlier_pick <- region_table$outlier_id

  # Count the 8 regions of interest (8)
nb <- length(outlier_pick) # count nb of regions 

  # Panel letters for each 8 regions of interest 
region_label_tibbles <- tibble(outlier_id = outlier_pick, label = letters[1:nb])

  # Set list of cool genes to find
cool_genes <-  c('kif16b', 'sox10', 
                 'casz1', 'hoxc8a','hoxc9','hoxc10a', 'hoxc13a',
                 'sws2abeta','sws2aalpha','sws2b','lws')

  # Set parameter of window buffer for plot
window_buffer <- 0.2*10^5

  # Set color parameter for plot line
gxp_clr <- c("#E69F00")

  # Set row count and measurement of each panels
row = 4
heights = c(6, 6)

  # Zoomed Manhattan plots 
p1 <- plot_regions(region_table, outlier_pick, region_label_tibbles, figure_path)
