
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


# Setting the current working directory (figure out a more robust way of doing this later)
setwd("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code")
source("./src/metrics.R")
source("./src/plotting.R")
source("./src/bootstrap.R")

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













