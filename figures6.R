
#' This is a file that would be used to outline the figures for chapter 6 in my thesis
#' that are made in R.
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")

################################################################################

# Figure - Custom PCA Plots

ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = "CLR")
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = "CLR")

# Plotting to check the correlation between the two datasets.
hex_colors <- blend_with_white(base_colours[conditions], timepts)

compsDNA <- pca.DNA$variates$X
compsRNA <- pca.RNA$variates$X


# Making my own PCA plots for DNA and RNA
df_compsDNA <- data.frame(compsDNA[,c(1,2)], "condition" = conditions, "colour" = hex_colors)
prop_expl_var <- pca.DNA$prop_expl_var$X
xlabel <- paste("PC1: ", round(100*prop_expl_var[1]), "% expl. var", sep = "")
ylabel <- paste("PC2: ", round(100*prop_expl_var[2]), "% expl. var", sep = "")

pca_plt1 <- df_compsDNA %>% ggplot(aes(x = PC1, y = PC2, shape = condition)) +
    geom_point(aes(color = colour), size = 3) +
    labs(title = "DNA", shape = "Condition",
         x = xlabel, y = ylabel) +
    scale_color_identity() +
    theme(legend.title = element_text("Condition")) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA), #"#FFFFFC"
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_textbox_simple(
              fill = "lightgray",
              linetype = 1,
              linewidth = 0.25,
              padding = margin(5, 5, 5, 5),
              halign = 0.5, valign = 0.5,
              face = "bold"
          ))

df_compsRNA <- data.frame(compsRNA[,c(1,2)], "condition" = conditions, "colour" = hex_colors)
prop_expl_var <- pca.RNA$prop_expl_var$X
xlabel <- paste("PC1: ", round(100*prop_expl_var[1]), "% expl. var", sep = "")
ylabel <- paste("PC2: ", round(100*prop_expl_var[2]), "% expl. var", sep = "")

pca_plt2 <- df_compsRNA %>% ggplot(aes(x = PC1, y = PC2, shape = condition)) +
    geom_point(aes(color = colour), size = 3) +
    labs(title = "RNA", shape = "Condition",
         x = xlabel, y = ylabel) +
    scale_color_identity() +
    theme(legend.title = element_text("Condition")) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_textbox_simple(
              fill = "lightgray",
              linetype = 1,
              linewidth = 0.25,
              padding = margin(5, 5, 5, 5),
              halign = 0.5, valign = 0.5,
              face = "bold"
          ))

plt <- (pca_plt1 + theme(legend.position = "none") | 
        pca_plt2 ) / wrapped_legend +
    plot_layout(heights = c(5,1.3)) +
    plot_annotation(theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))

ggsave("PCA.png", plot = plt, path = "./figures/",
       width = 24, height = 14, units = "cm")



################################################################################

# Plotting the Bray-Curtis, PCA Bray-Curtis (λ = 0 and 1)

# Original Bray-Curtis
DNABray <- sampledist(DNA.1, Bray_Curtis)
title <- expression(bold("DNA"))

MDSstructure <- mdsplotfunc(DNABray, conditions, timepts, title)
MDSstructure$stress

MDS_BC_DNAplt <- MDSstructure$mdsplot
MDS_BC_DNAplt <- MDS_BC_DNAplt + theme(axis.title.x = element_blank(), 
                                 legend.position = "none") +
    ylab(expression(atop(bold("Bray-Curtis"), "D2"))) 
    


RNABray <- sampledist(RNA.1, Bray_Curtis)
title <- expression(bold("RNA"))

MDSstructure <- mdsplotfunc(RNABray, conditions, timepts, title)
MDSstructure$stress

MDS_BC_RNAplt <- MDSstructure$mdsplot
MDS_BC_RNAplt <- MDS_BC_RNAplt + theme(axis.title.x = element_blank(), 
                                 axis.title.y = element_blank(),
                                 legend.position = "none")



# DNA: PCA Bary-Curtis with and without reactor time-scale
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
DNAcomps <- pca.DNA$variates$X

DNA_PCABray0 <- sampledist(DNAcomps, PCA_Bray_Curtis)
title <- expression(bold("DNA"))

MDSstructure <- mdsplotfunc(DNA_PCABray0, conditions, timepts, title)
MDSstructure$stress

MDS_PCABC0_DNAplt <- MDSstructure$mdsplot
MDS_PCABC0_DNAplt <- MDS_PCABC0_DNAplt + theme(axis.title.x = element_blank(), 
                                               legend.position = "none",
                                               plot.title = element_blank()) +
    ylab(expression(atop(bold("PCA Bray-Curtis (λ = 0)"), "D2")))


lambda <- 1
DNA_PCABray1 <- penaltyfunc(df = DNAcomps,
                            metric = PCA_Bray_Curtis,
                            reactor = reactors,
                            timepts = timepts,
                            kernelfunc = linearkernel,
                            lambda = lambda)
title <- expression(bold("DNA"))

MDSstructure <- mdsplotfunc(DNA_PCABray1, conditions, timepts, title)
MDSstructure$stress

MDS_PCABC1_DNAplt <- MDSstructure$mdsplot
MDS_PCABC1_DNAplt <- MDS_PCABC1_DNAplt + theme(legend.position = "none",
                                               plot.title = element_blank()) +
    ylab(expression(atop(bold("PCA Bray-Curtis (λ = 1)"), "D2")))



# RNA: PCA Bary-Curtis with and without reactor time-scale
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
RNAcomps <- pca.RNA$variates$X

RNA_PCABray0 <- sampledist(RNAcomps, PCA_Bray_Curtis)
title <- expression(bold("RNA"))

MDSstructure <- mdsplotfunc(RNA_PCABray0, conditions, timepts, title)
MDSstructure$stress

MDS_PCABC0_RNAplt <- MDSstructure$mdsplot
MDS_PCABC0_RNAplt <- MDS_PCABC0_RNAplt + theme(axis.title.x = element_blank(), 
                                               axis.title.y = element_blank(),
                                               plot.title = element_blank()) 


lambda <- 1
RNA_PCABray1 <- penaltyfunc(df = RNAcomps,
                            metric = PCA_Bray_Curtis,
                            reactor = reactors,
                            timepts = timepts,
                            kernelfunc = linearkernel,
                            lambda = lambda)
title <- expression(bold("RNA"))

MDSstructure <- mdsplotfunc(RNA_PCABray1, conditions, timepts, title)
MDSstructure$stress

MDS_PCABC1_RNAplt <- MDSstructure$mdsplot
MDS_PCABC1_RNAplt <- MDS_PCABC1_RNAplt + theme(legend.position = "none",
                                               axis.title.y = element_blank(),
                                               plot.title = element_blank()) 


plt <- (MDS_BC_DNAplt | MDS_BC_RNAplt) /
    (MDS_PCABC0_DNAplt | MDS_PCABC0_RNAplt) /
    (MDS_PCABC1_DNAplt | MDS_PCABC1_RNAplt) / 
    wrapped_legend + 
    plot_layout(heights = c(4,4,4,1.5))

ggsave("singleBray.png", plot = plt, path = "./figures/",
       width = 24, height = 28, units = "cm")

################################################################################

# Numerics on the MDS plots from the PCA Bray-Curtis in last chunk

# DNA
DNA_BrayNum <- MDSnumeric(df = DNABray,
                      func = Bray_Curtis,
                      reactor = reactors,
                      condition = conditions,
                      timepts = timepts,
                      distdf = TRUE)
DNA_BrayNum


DNA_PCABrayNum0 <- MDSnumeric(df = DNA_PCABray0,
                          func = PCA_Bray_Curtis,
                          reactor = reactors,
                          condition = conditions,
                          timepts = timepts,
                          distdf = TRUE)
DNA_PCABrayNum0


DNA_PCABrayNum1 <- MDSnumeric(df = DNA_PCABray1,
                          func = PCA_Bray_Curtis,
                          reactor = reactors,
                          condition = conditions,
                          timepts = timepts,
                          distdf = TRUE)
DNA_PCABrayNum1


# RNA
RNA_BrayNum <- MDSnumeric(df = RNABray,
                          func = Bray_Curtis,
                          reactor = reactors,
                          condition = conditions,
                          timepts = timepts,
                          distdf = TRUE)
RNA_BrayNum


RNA_PCABrayNum0 <- MDSnumeric(df = RNA_PCABray0,
                              func = PCA_Bray_Curtis,
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
RNA_PCABrayNum0


RNA_PCABrayNum1 <- MDSnumeric(df = RNA_PCABray1,
                              func = PCA_Bray_Curtis,
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
RNA_PCABrayNum1

################################################################################
# Figure - Custom PLS plot

ncomps <- 4
X1 <- logratio.transfo(DNA.1 + 0.001, logratio = "CLR")
X2 <- logratio.transfo(RNA.1 + 0.001, logratio = "CLR")

pls.result <- pls(X=X1, Y=X2, ncomp = ncomps, mode = "canonical")

# Plotting to check the correlation between the two datasets.
hex_colors <- blend_with_white(base_colours[conditions], timepts)


compsDNA <- pls.result$variates$X
compsRNA <- pls.result$variates$Y
rownames(compsRNA) <- rownames(RNA.1)


# Making my own PLS plots for DNA and RNA
df_compsDNA <- data.frame(compsDNA[,c(1,2)], "condition" = conditions, "colour" = hex_colors)
prop_expl_var <- pls.result$prop_expl_var$X
xlabel <- paste("Variate 1: ", round(100*prop_expl_var[1]), "% expl. var", sep = "")
ylabel <- paste("Variate 2: ", round(100*prop_expl_var[2]), "% expl. var", sep = "")

pls_plt1 <- df_compsDNA %>% ggplot(aes(x = comp1, y = comp2, shape = condition)) +
    geom_point(aes(color = colour), size = 3) +
    labs(title = "X = DNA", shape = "Condition",
         x = xlabel, y = ylabel) +
    scale_color_identity() +
    theme(legend.title = element_text("Condition")) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_textbox_simple(
              fill = "lightgray",
              linetype = 1,
              linewidth = 0.25,
              padding = margin(5, 5, 5, 5),
              halign = 0.5, valign = 0.5,
              face = "bold"
          ))

df_compsRNA <- data.frame(compsRNA[,c(1,2)], "condition" = conditions, "colour" = hex_colors)
prop_expl_var <- pls.result$prop_expl_var$Y
xlabel <- paste("Variate 1: ", round(100*prop_expl_var[1]), "% expl. var", sep = "")
ylabel <- paste("Variate 2: ", round(100*prop_expl_var[2]), "% expl. var", sep = "")

pls_plt2 <- df_compsRNA %>% ggplot(aes(x = comp1, y = comp2, shape = condition)) +
    geom_point(aes(color = colour), size = 3) +
    labs(title = "Y = RNA", shape = "Condition",
         x = xlabel, y = ylabel) +
    scale_color_identity() +
    theme(legend.title = element_text("Condition")) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_textbox_simple(
              fill = "lightgray",
              linetype = 1,
              linewidth = 0.25,
              padding = margin(5, 5, 5, 5),
              halign = 0.5, valign = 0.5,
              face = "bold"
          ))

plt <- (pls_plt1 + theme(legend.position = "none") | 
        pls_plt2) / wrapped_legend +
    plot_layout(heights = c(5,1.3)) +
    plot_annotation(theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))
plt

ggsave("PLS.png", plot = plt, path = "./figures/",
       width = 24, height = 14, units = "cm")



# Using the PLS results, calculate the correlation values between the DNA and RNA components
# Just make sure you are using the PLS components and not the PCA ones.
diag(cor(compsDNA, compsRNA))

################################################################################

# comparing the PLS and PCA Bray-Curtis (no time-scale)

# PCA
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = "CLR")
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = "CLR")

PCAcompsDNA <- pca.DNA$variates$X
PCAcompsRNA <- pca.RNA$variates$X

#PLS
ncomps <- 4
X1 <- logratio.transfo(DNA.1 + 0.001, logratio = "CLR")
X2 <- logratio.transfo(RNA.1 + 0.001, logratio = "CLR")

pls.result <- pls(X=X1, Y=X2, ncomp = ncomps, mode = "canonical")
PLScompsDNA <- pls.result$variates$X
PLScompsRNA <- pls.result$variates$Y
rownames(PLScompsRNA) <- rownames(RNA.1)


#NDMS plots
DNA_PCA <- sampledist(PCAcompsDNA, PCA_Bray_Curtis)
title <- expression(bold("DNA"))

MDSstructure <- mdsplotfunc(DNA_PCA, conditions, timepts, title)
MDSstructure$stress

MDS_DNA_PCA <- MDSstructure$mdsplot
MDS_DNA_PCA <- MDS_DNA_PCA + theme(axis.title.x = element_blank(), 
                                               legend.position = "none") +
    ylab(expression(atop(bold("PCA Bray-Curtis"), "D2")))

##
RNA_PCA <- sampledist(PCAcompsRNA, PCA_Bray_Curtis)
title <- expression(bold("RNA"))

MDSstructure <- mdsplotfunc(RNA_PCA, conditions, timepts, title)
MDSstructure$stress

MDS_RNA_PCA <- MDSstructure$mdsplot
MDS_RNA_PCA <- MDS_RNA_PCA + theme(axis.title.x = element_blank(), 
                                   axis.title.y = element_blank(),
                                   legend.position = "none") 


##
DNA_PLS <- sampledist(PLScompsDNA, PCA_Bray_Curtis)
title <- expression(bold("DNA"))

MDSstructure <- mdsplotfunc(DNA_PLS, conditions, timepts, title)
MDSstructure$stress

MDS_DNA_PLS <- MDSstructure$mdsplot
MDS_DNA_PLS <- MDS_DNA_PLS + theme(plot.title = element_blank(),
                                   legend.position = "none") +
    ylab(expression(atop(bold("PLS Bray-Curtis"), "D2")))

##
RNA_PLS <- sampledist(PLScompsRNA, PCA_Bray_Curtis)
title <- expression(bold("RNA"))

MDSstructure <- mdsplotfunc(RNA_PLS, conditions, timepts, title)
MDSstructure$stress

MDS_RNA_PLS <- MDSstructure$mdsplot
MDS_RNA_PLS <- MDS_RNA_PLS + theme(axis.title.y = element_blank(),
                                   plot.title = element_blank(),
                                   legend.position = c(1.15,1.05)) 


## 
plt <- (MDS_DNA_PCA | MDS_RNA_PCA) /
    (MDS_DNA_PLS | MDS_RNA_PLS) /
    wrapped_legend + 
    plot_layout(heights = c(4,4,1.5)) + 
    plot_annotation(theme = theme(plot.margin = unit(c(0.5, 3, 0, 0), "cm")))

ggsave("PCAvsPLS.png", plot = plt, path = "./figures/",
       width = 24.1, height = 21, units = "cm")


################################################################################

# NMDS on DNA, RNA and Integrated (for PLS Bray-Curtis with and without time-scale)

#PLS
ncomps <- 4
X1 <- logratio.transfo(DNA.1 + 0.001, logratio = "CLR")
X2 <- logratio.transfo(RNA.1 + 0.001, logratio = "CLR")

pls.result <- pls(X=X1, Y=X2, ncomp = ncomps, mode = "canonical")
PLScompsDNA <- pls.result$variates$X
PLScompsRNA <- pls.result$variates$Y
rownames(PLScompsRNA) <- rownames(RNA.1)


# DNA
DNA_PLSBray0 <- sampledist(PLScompsDNA, PCA_Bray_Curtis)
title <- expression(bold("DNA"))

MDSstructure <- mdsplotfunc(DNA_PLSBray0, conditions, timepts, title)
MDSstructure$stress

MDS_PLSBC0_DNAplt <- MDSstructure$mdsplot
MDS_PLSBC0_DNAplt <- MDS_PLSBC0_DNAplt + theme(axis.title.x = element_blank(), 
                                               legend.position = "none") +
    ylab(expression(atop(bold("PLS Bray-Curtis (λ = 0)"), "D2")))


lambda <- 1
DNA_PLSBray1 <- penaltyfunc(df = PLScompsDNA,
                            metric = PCA_Bray_Curtis,
                            reactor = reactors,
                            timepts = timepts,
                            kernelfunc = linearkernel,
                            lambda = lambda)
title <- expression(bold("DNA"))

MDSstructure <- mdsplotfunc(DNA_PLSBray1, conditions, timepts, title)
MDSstructure$stress

MDS_PLSBC1_DNAplt <- MDSstructure$mdsplot
MDS_PLSBC1_DNAplt <- MDS_PLSBC1_DNAplt + theme(legend.position = "none",
                                               plot.title = element_blank()) +
    ylab(expression(atop(bold("PLS Bray-Curtis (λ = 1)"), "D2")))

# RNA
RNA_PLSBray0 <- sampledist(PLScompsRNA, PCA_Bray_Curtis)
title <- expression(bold("RNA"))

MDSstructure <- mdsplotfunc(RNA_PLSBray0, conditions, timepts, title)
MDSstructure$stress

MDS_PLSBC0_RNAplt <- MDSstructure$mdsplot
MDS_PLSBC0_RNAplt <- MDS_PLSBC0_RNAplt + theme(axis.title.x = element_blank(), 
                                               legend.position = "none") 


lambda <- 1
RNA_PLSBray1 <- penaltyfunc(df = PLScompsRNA,
                            metric = PCA_Bray_Curtis,
                            reactor = reactors,
                            timepts = timepts,
                            kernelfunc = linearkernel,
                            lambda = lambda)
title <- expression(bold("RNA"))

MDSstructure <- mdsplotfunc(RNA_PLSBray1, conditions, timepts, title)
MDSstructure$stress

MDS_PLSBC1_RNAplt <- MDSstructure$mdsplot
MDS_PLSBC1_RNAplt <- MDS_PLSBC1_RNAplt + theme(legend.position = "none",
                                               plot.title = element_blank()) 

# Integrated
INTcomps <- (PLScompsDNA + PLScompsRNA) / 2

INT_PLSBray0 <- sampledist(INTcomps, PCA_Bray_Curtis)
title <- expression(bold("Integrated"))

MDSstructure <- mdsplotfunc(INT_PLSBray0, conditions, timepts, title)
MDSstructure$stress

MDS_PLSBC0_INTplt <- MDSstructure$mdsplot
MDS_PLSBC0_INTplt <- MDS_PLSBC0_INTplt + theme(axis.title.x = element_blank(), 
                                               legend.position = "none") 


lambda <- 1
INT_PLSBray1 <- penaltyfunc(df = INTcomps,
                            metric = PCA_Bray_Curtis,
                            reactor = reactors,
                            timepts = timepts,
                            kernelfunc = linearkernel,
                            lambda = lambda)
title <- expression(bold("Integrated"))

MDSstructure <- mdsplotfunc(INT_PLSBray1, conditions, timepts, title)
MDSstructure$stress

MDS_PLSBC1_INTplt <- MDSstructure$mdsplot
MDS_PLSBC1_INTplt <- MDS_PLSBC1_INTplt + theme(legend.position = c(1.2,1.05),
                                               plot.title = element_blank()) 



# Combine the plots
plt <- (MDS_PLSBC0_DNAplt | MDS_PLSBC0_RNAplt | MDS_PLSBC0_INTplt) /
    (MDS_PLSBC1_DNAplt | MDS_PLSBC1_RNAplt | MDS_PLSBC1_INTplt) +
    wrapped_legend +
    plot_layout(heights = c(4,4,1.5)) +
    plot_annotation(theme = theme(plot.margin = unit(c(0.5, 3, 0, 0), "cm")))

ggsave("multiBray.png", plot = plt, path = "./figures/",
       width = 32, height = 21, units = "cm")



################################################################################

# Numerics on the MDS plots from the PLS comparison in the last chunk

# DNA
DNA_PLSBrayNum0 <- MDSnumeric(df = DNA_PLSBray0,
                          func = Bray_Curtis,
                          reactor = reactors,
                          condition = conditions,
                          timepts = timepts,
                          distdf = TRUE)
DNA_PLSBrayNum0


DNA_PLSBrayNum1 <- MDSnumeric(df = DNA_PLSBray1,
                              func = PCA_Bray_Curtis,
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
DNA_PLSBrayNum1


# RNA
RNA_PLSBrayNum0 <- MDSnumeric(df = RNA_PLSBray0,
                              func = PCA_Bray_Curtis,
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
RNA_PLSBrayNum0



RNA_PLSBray1Num <- MDSnumeric(df = RNA_PLSBray1,
                          func = Bray_Curtis,
                          reactor = reactors,
                          condition = conditions,
                          timepts = timepts,
                          distdf = TRUE)
RNA_PLSBray1Num


# Integrated
INT_PLSBrayNum0 <- MDSnumeric(df = INT_PLSBray0,
                              func = PCA_Bray_Curtis,
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
INT_PLSBrayNum0


INT_PLSBrayNum1 <- MDSnumeric(df = INT_PLSBray1,
                              func = PCA_Bray_Curtis,
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
INT_PLSBrayNum1



























