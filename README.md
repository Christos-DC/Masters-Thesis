# Thesis Title - Ecological Distance Metrics for Temporal Omics Data in Microbiome Studies

**Supervisors**: Kim-Anh L&ecirc; Cao, Saritha Kodikara

**Date**: 15/05/2026

**Schools**: Melbourne Integrative Genomics (MIG), School of Mathematics and Statistics

**Institute**: The University of Melbourne


## Abstract
Microbiome studies explain how bacterial species change when exposed to different environmental conditions. Changes in bacterial composition are captured by using their omics features like DNA and RNA to quantify information like species abundance in the community. A prevalent interest in these studies is comparing species diversity across different microbial communities, which are calculated by using ecological distance metrics.

Ecological distance metrics are tools that measure diversity for different bacteria sample groups. Many distance metrics have been proposed by microbiologists to capture similarity between samples. However, these metrics assume that samples are independent
from each other and are based on one omics dataset, which does not hold when samples are recorded over time or when multiple omics datasets are present.

The focus of this work was on developing a new metric that explains the diversity evolution of species within samples and can be used for multi-omics integration. Temporal correlations are first captured using ordination methods such as Principal Component
Analysis and Partial Least Squares, and are then incorporated into this metric for both single and multi-omics frameworks.

Motivated by an anaerobic digestion study, this work analyses how anaerobic digestion under salt exposure affects species abundance over a 14-week time frame. Abundance was recorded using the omics features DNA and RNA. The new metric was first used on each of these features, then later implemented in an integrated multi-omics framework to capture temporal changes over multiple biological layers. The results show that this approach provided a more complete view of biological processes and detected coordinated responses better compared to existing distance metrics.


## File Structure
The code was developed in R under version 4.5.3. Ensure that your version of R is at least 4.5.3 or above for the code to run properly.


