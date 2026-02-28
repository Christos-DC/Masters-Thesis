
# Loading the Packages needed for this project

# Way to install the mixOmics and phlyoseq packages if not installed already
if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install("mixOmics")
    BiocManager::install("phyloseq")
}

# mixOmics
library(mixOmics)
library(phyloseq)

# Parallel execution packages
library(foreach)
library(doParallel)

# Data manipulation and visualisations
library(ggplot2)
library(dplyr)
library(readr)
library(smacof) 
library(scales)
library(patchwork)

# Installing Packages for Heatmaps via github
# 
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)


# Gather the source R codes
setwd("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code")
source("src/metrics.R", local = TRUE)
source("src/plotting.R", local = TRUE)
source("src/bootstrap.R", local = TRUE)
source("src/kernels.R", local = TRUE)

############################### Data Loading ###################################
DNA <- readRDS("./zymo_DNA.rds")
RNA <- readRDS("./zymo_RNA.rds")

DNA.1 <- phyloseq(otu_table(DNA))
DNA.2 <- phyloseq(tax_table(DNA))
DNA.3 <- phyloseq(sample_data(DNA))

DNA.1 <- data.frame(DNA.1)
DNA.2 <- data.frame(DNA.2)
DNA.3 <- data.frame(DNA.3)

RNA.1 <- phyloseq(otu_table(RNA))
RNA.2 <- phyloseq(tax_table(RNA))
RNA.3 <- phyloseq(sample_data(RNA))

RNA.1 <- data.frame(RNA.1)
RNA.2 <- data.frame(RNA.2)
RNA.3 <- data.frame(RNA.3)

phylonames <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(DNA.2) <- phylonames
colnames(RNA.2) <- phylonames

############################## Pre-filtering ###################################

low.count.removal <- function(data, percent = 0.01){
    #' Inputs:
    #'      data (data.frame): OTU count df of size n (sample) x p (OTU)
    #'      percent (float): Positive integer for threshold in rejecting very 
    #'                       small OTU counts.
    #' 
    #' Returns: A list containing two items
    #'      data.filter (data.frame): Filtered version of the data frame.
    #'      keep.otu (vector): A binary vector outlining which OTU's to keep.
    
    keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
    data.filter = data[, keep.otu]
    return(list(data.filter = data.filter, keep.otu = keep.otu))
}

DNA.1 <- t(DNA.1)
DNA.1 <- low.count.removal(DNA.1)$data.filter

RNA.1 <- t(RNA.1)
RNA.1 <- low.count.removal(RNA.1)$data.filter


########################## Legend Plot for MDS #################################

# Plotting the legend plot for Control, Inhibition 1 and 2.
mintime <- min(DNA.3$week_from_salt1)
maxtime <- max(DNA.3$week_from_salt1)
n <- 100
weeks <- round(seq(mintime, maxtime, length.out = n), 4)

# Gather the blended colours over time.
legend_df <- expand.grid(sort(names(base_colours)), weeks)
colnames(legend_df) <- c("condition", "week")
legend_df$color <- blend_with_white(base_colours[legend_df$condition], as.numeric(legend_df$week)) 

# Legend plot (compact, horizontal)
legend_plot <- ggplot(legend_df, aes(x = week, y = condition, fill = color)) +
    geom_tile(height = 0.6) +
    scale_fill_identity() +
    theme_minimal(base_size = 7) +
    labs(x = "Week", y = NULL) +
    scale_x_continuous(breaks = seq(-3,11, by=1)) +
    theme(
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", color = NA)
    )

# Control the dimension of the plot
wrapped_legend <- wrap_elements(legend_plot) +
    plot_layout(widths = unit(0.5, "npc"))









