

# Figures for Oliver to see
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")

################################################################################
# My PCA plots
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = "CLR")
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = "CLR")

# Plotting to check the correlation between the two datasets.
hex_colors <- blend_with_white(base_colours[conditions], timepts)

compsDNA <- pca.DNA$variates$X
compsRNA <- pca.RNA$variates$X


# Making my own PCA plots for DNA and RNA
df_compsDNA <- data.frame(compsDNA[,-c(3,4)], "condition" = conditions, "colour" = hex_colors)
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
              linewidth = 0.2,
              padding = margin(5, 5, 5, 5),
              halign = 0.5, valign = 0.5
          ))

df_compsRNA <- data.frame(compsRNA[,-c(3,4)], "condition" = conditions, "colour" = hex_colors)
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
              linewidth = 0.2,
              padding = margin(5, 5, 5, 5),
              halign = 0.5, valign = 0.5
          ))

(pca_plt1 + theme(legend.position = "none") | 
        pca_plt2 ) / wrapped_legend +
    plot_layout(heights = c(5,1.3)) +
    plot_annotation(title = "PCA plots for DNA and RNA", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))






# DM Bray Curtis on DNA and RNA plots 
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.DNA$variates$X[, 1:ncomps]

DNA_PCABray <- sampledist(comps, PCA_Bray_Curtis)
title <- "DNA"

MDSstructure <- mdsplotfunc(DNA_PCABray, conditions, timepts, title)
MDSstructure$stress

DNA_MDSplt <- MDSstructure$mdsplot


# PCA Bray-Curtis (RNA)
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.RNA$variates$X[, 1:ncomps]

RNA_PCABray <- sampledist(comps, PCA_Bray_Curtis)
title <- "RNA"

MDSstructure <- mdsplotfunc(RNA_PCABray, conditions, timepts, title)
MDSstructure$stress

RNA_MDSplt <- MDSstructure$mdsplot

plt <- (DNA_MDSplt + theme(legend.position = "none")| RNA_MDSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5,1.5)) +
    plot_annotation(title = "NMDS using 4 PCA components for DNA and RNA datasets", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))
plt


################################################################################
# DM Bray Curtis on DNA and RNA with time adjustment plots

lambda <- 1
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.DNA$variates$X[, 1:ncomps]

text <- "DNA"
DNAadjust <- penaltyfunc(comps, PCA_Bray_Curtis, reactors, timepts, linearkernel, lambda)
mdsresult <- mdsplotfunc(DNAadjust,conditions, timepts, text)

mdsresult$stress
DNA_MDSplt <- mdsresult$mdsplot


pca.RNA <- pca(RNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.RNA$variates$X[, 1:ncomps]

text <- "RNA"
RNAadjust <- penaltyfunc(comps, PCA_Bray_Curtis, reactors, timepts, linearkernel, lambda)
mdsresult <- mdsplotfunc(RNAadjust,conditions, timepts, text)

mdsresult$stress
RNA_MDSplt <- mdsresult$mdsplot




plt <- (DNA_MDSplt + theme(legend.position = "none")| RNA_MDSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5,1.5)) +
    plot_annotation(title = "NMDS using 4 PCA components for DNA and RNA datasets:\n Linear Kernel (λ = 1)", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))
plt


################################################################################
# N-integration step on DNA and RNA

# Making the PLS plot
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
df_compsDNA <- data.frame(compsDNA[,-c(3,4)], "condition" = conditions, "colour" = hex_colors)

pls_plt1 <- df_compsDNA %>% ggplot(aes(x = comp1, y = comp2, shape = condition)) +
    geom_point(aes(color = colour), size = 3) +
    labs(title = "X = DNA", shape = "Condition",
         x = "variate 1", y = "variate 2") +
    scale_color_identity() +
    theme(legend.title = element_text("Condition")) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "#FFFFFC", color = NA),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(size = 12, hjust = 0.5))

df_compsRNA <- data.frame(compsRNA[,-c(3,4)], "condition" = conditions, "colour" = hex_colors)

pls_plt2 <- df_compsRNA %>% ggplot(aes(x = comp1, y = comp2, shape = condition)) +
    geom_point(aes(color = colour), size = 3) +
    labs(title = "Y = RNA", shape = "Condition",
         x = "variate 1", y = "variate 2") +
    scale_color_identity() +
    theme(legend.title = element_text("Condition")) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "#FFFFFC", color = NA),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(size = 12, hjust = 0.5))

(pls_plt1 + theme(legend.position = "none") | 
        pls_plt2 + theme(axis.title.y = element_blank())) / wrapped_legend +
    plot_layout(heights = c(5,1.3)) +
    plot_annotation(title = "PLS plots for DNA and RNA", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))






# No Time adjustment
## For the DNA component
DNA_PLSBray <- sampledist(compsDNA, PCA_Bray_Curtis)
title <- "DNA"

MDSstructure <- mdsplotfunc(DNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

DNA_MDSplt <- MDSstructure$mdsplot

## For the RNA component
RNA_PLSBray <- sampledist(compsRNA, PCA_Bray_Curtis)
title <- "RNA"

MDSstructure <- mdsplotfunc(RNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

RNA_MDSplt <- MDSstructure$mdsplot



plt <- (DNA_MDSplt + theme(legend.position = "none")| RNA_MDSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5,1.5)) +
    plot_annotation(title = "NMDS using 4 PLS components for DNA and RNA datasets before Integration", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))
plt



## Averaging the components
compsCombined <- (compsDNA + compsRNA)/2


Comb_PLSBray <- sampledist(compsCombined, PCA_Bray_Curtis)
title <- "Integrated"

MDSstructure <- mdsplotfunc(Comb_PLSBray, conditions, timepts, title)
MDSstructure$stress

Int_MDSplt <- MDSstructure$mdsplot
Int_MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3)) +
    plot_annotation(title = "NMDS using 4 PLS components for DNA and RNA datasets after Integration", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))


# Showing the three plots side-by-side
(DNA_MDSplt + theme(legend.position = "none")| 
        RNA_MDSplt + theme(axis.title.y = element_blank(), legend.position = "none") | 
        Int_MDSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5,1.5)) +
    plot_annotation(title = "NMDS using 4 PLS components for DNA, RNA and Integrated Plots", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))







# Time adjustment
lambda <- 1
DNA_PLSBray <- penaltyfunc(df = compsDNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
title <- "DNA"

MDSstructure <- mdsplotfunc(DNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

DNA_MDSplt <- MDSstructure$mdsplot



# For the RNA component
RNA_PLSBray <- penaltyfunc(df = compsRNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
title <- "RNA"

MDSstructure <- mdsplotfunc(RNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

RNA_MDSplt <- MDSstructure$mdsplot



Comb_PLSBray <- penaltyfunc(df = compsCombined,
                            metric = PCA_Bray_Curtis,
                            reactor = reactors,
                            timepts = timepts,
                            kernelfunc = linearkernel,
                            lambda = lambda)
title <- "Integrated"

MDSstructure <- mdsplotfunc(Comb_PLSBray, conditions, timepts, title)
MDSstructure$stress

Int_MDSplt <- MDSstructure$mdsplot

(DNA_MDSplt + theme(legend.position = "none")| 
        RNA_MDSplt + theme(axis.title.y = element_blank(), legend.position = "none") | 
        Int_MDSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5,1.5)) +
    plot_annotation(title = "NMDS using 4 PLS components for DNA, RNA and Integrated Plots:\n Linear Kernel (λ = 1)", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))


################################################################################

#### Trying to set up the Bray Curtis, PCA Bray-Curtis (lam = 0 and lam = 1) ####

# Normal Bray-Curtis
BC_DNA <- sampledist(DNA.1, Bray_Curtis)
title <- expression(bold("Bray-Curtis"))

MDSstructure <- mdsplotfunc(BC_DNA, conditions, timepts, title)
MDSstructure$stress

BC_DNA_MDSplt <- MDSstructure$mdsplot
BC_DNA_MDSplt <- BC_DNA_MDSplt + 
    ylab(expression(atop(bold("DNA"), "D2"))) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          plot.title = element_text(size = 14))
BC_DNA_MDSplt

#
BC_RNA <- sampledist(RNA.1, Bray_Curtis)
title <- expression(bold("Bray-Curtis"))

MDSstructure <- mdsplotfunc(BC_RNA, conditions, timepts, title)
MDSstructure$stress

BC_RNA_MDSplt <- MDSstructure$mdsplot
BC_RNA_MDSplt <- BC_RNA_MDSplt + 
    ylab(expression(atop(bold("RNA"), "D2"))) +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          plot.title = element_blank())
BC_RNA_MDSplt


# DM Bray Curtis on DNA and RNA plots 
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.DNA$variates$X[, 1:ncomps]

DNA_PCABray <- sampledist(comps, PCA_Bray_Curtis)
title <- expression(bold("PCA Bray-Curtis"))

MDSstructure <- mdsplotfunc(DNA_PCABray, conditions, timepts, title)
MDSstructure$stress

PCA_DNA_MDSplt0 <- MDSstructure$mdsplot
PCA_DNA_MDSplt0 <- PCA_DNA_MDSplt0 + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 14))
PCA_DNA_MDSplt0


# PCA Bray-Curtis (RNA)
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.RNA$variates$X[, 1:ncomps]

RNA_PCABray <- sampledist(comps, PCA_Bray_Curtis)
title <- "RNA"

MDSstructure <- mdsplotfunc(RNA_PCABray, conditions, timepts, title)
MDSstructure$stress

PCA_RNA_MDSplt0 <- MDSstructure$mdsplot
PCA_RNA_MDSplt0 <- PCA_RNA_MDSplt0 + 
    theme(legend.position = "none",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_blank(),
          plot.title = element_blank())
PCA_RNA_MDSplt0



#DM Bray Curtis on DNA and RNA with time adjustment λ=1
lambda <- 1
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.DNA$variates$X[, 1:ncomps]

DNA_PCABray <- penaltyfunc(df = comps,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
title <- expression(bold("PCA Bray-Curtis with Time-Scale (λ=1)"))

MDSstructure <- mdsplotfunc(DNA_PCABray, conditions, timepts, title)
MDSstructure$stress

PCA_DNA_MDSplt1 <- MDSstructure$mdsplot
PCA_DNA_MDSplt1 <- PCA_DNA_MDSplt1 + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 14))
PCA_DNA_MDSplt1

#
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.RNA$variates$X[, 1:ncomps]

RNA_PCABray <- penaltyfunc(df = comps,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
title <- expression(bold("PCA Bray-Curtis with Time-Scale (λ=1)"))

MDSstructure <- mdsplotfunc(RNA_PCABray, conditions, timepts, title)
MDSstructure$stress

PCA_RNA_MDSplt1 <- MDSstructure$mdsplot
PCA_RNA_MDSplt1 <- PCA_RNA_MDSplt1 + 
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_blank(),
          plot.title = element_blank(),
          legend.position = c(1.15,1.1))
PCA_RNA_MDSplt1



(BC_DNA_MDSplt | PCA_DNA_MDSplt0 | PCA_DNA_MDSplt1) /
    (BC_RNA_MDSplt | PCA_RNA_MDSplt0 | PCA_RNA_MDSplt1) /
    wrapped_legend +
    plot_layout(heights = c(4,4,2)) +
    plot_annotation(title = "MDS Plots: Comparing Normal Bray-Curtis, PCA Bray-Curtis With and Without Time Scale\n", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5),
                                  plot.margin = unit(c(0.5, 3, 0, 0), "cm")))


##################################################################################
# Get PLS components
ncomps <- 4
X1 <- logratio.transfo(DNA.1 + 0.001, logratio = "CLR")
X2 <- logratio.transfo(RNA.1 + 0.001, logratio = "CLR")

pls.result <- pls(X=X1, Y=X2, ncomp = ncomps, mode = "canonical")

PLScompsDNA <- pls.result$variates$X
PLScompsRNA <- pls.result$variates$Y

rownames(PLScompsRNA) <- rownames(RNA.1)


# DM Bray Curtis on DNA and RNA plots 
DNA_PLSBray <- sampledist(PLScompsDNA, PCA_Bray_Curtis)
title <- expression(bold("DNA"))

MDSstructure <- mdsplotfunc(DNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

PLS_DNA_MDSplt0 <- MDSstructure$mdsplot
PLS_DNA_MDSplt0 <- PLS_DNA_MDSplt0 + 
    ylab(expression(atop(bold("PLS Bray-Curtis"), "D2"))) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          plot.title = element_text(size = 12))
PLS_DNA_MDSplt0


# PCA Bray-Curtis (RNA)
RNA_PLSBray <- sampledist(PLScompsRNA, PCA_Bray_Curtis)
title <- expression(bold("RNA"))

MDSstructure <- mdsplotfunc(RNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

PLS_RNA_MDSplt0 <- MDSstructure$mdsplot
PLS_RNA_MDSplt0 <- PLS_RNA_MDSplt0 + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 12))
PLS_RNA_MDSplt0



#DM Bray Curtis on DNA and RNA with time adjustment λ=1
lambda <- 1
DNA_PLSBray <- penaltyfunc(df = PLScompsDNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
title <- expression(bold("PLS Bray-Curtis with Time-Scale (λ=1)"))

MDSstructure <- mdsplotfunc(DNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

PLS_DNA_MDSplt1 <- MDSstructure$mdsplot
PLS_DNA_MDSplt1 <- PLS_DNA_MDSplt1 + 
    ylab(expression(atop(bold("PLS Bray-Curtis with Time-Scale (λ=1)"), "D2"))) +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          plot.title = element_blank())
PLS_DNA_MDSplt1

#
RNA_PLSBray <- penaltyfunc(df = PLScompsRNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
title <- expression(bold("PLS Bray-Curtis with Time-Scale (λ=1)"))

MDSstructure <- mdsplotfunc(RNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

PLS_RNA_MDSplt1 <- MDSstructure$mdsplot
PLS_RNA_MDSplt1 <- PLS_RNA_MDSplt1 + 
    theme(legend.position = "none",
          axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(),
          plot.title = element_blank())
PLS_RNA_MDSplt1



# Doing the integration step

## Combining the components together
compsCombined <- (PLScompsDNA + PLScompsRNA)/2
INTG_PLSBray <- sampledist(compsCombined, PCA_Bray_Curtis)

title <- expression(bold("Integrated"))

MDSstructure <- mdsplotfunc(INTG_PLSBray, conditions, timepts, title)
MDSstructure$stress

PLS_INTG_MDSplt0 <- MDSstructure$mdsplot
PLS_INTG_MDSplt0 <- PLS_INTG_MDSplt0 + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 12))
PLS_INTG_MDSplt0


## With lambda = 1
lambda <- 1
compsCombined <- (PLScompsDNA + PLScompsRNA)/2
INTG_PLSBray <- penaltyfunc(df = compsCombined,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)

title <- expression(bold("Integrated"))

MDSstructure <- mdsplotfunc(INTG_PLSBray, conditions, timepts, title)
MDSstructure$stress

PLS_INTG_MDSplt1 <- MDSstructure$mdsplot
PLS_INTG_MDSplt1 <- PLS_INTG_MDSplt1 + 
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(),
          plot.title = element_blank(),
          legend.position = c(1.15,1.1))
PLS_INTG_MDSplt1


(PLS_DNA_MDSplt0 | PLS_RNA_MDSplt0 | PLS_INTG_MDSplt0) /
    (PLS_DNA_MDSplt1 | PLS_RNA_MDSplt1 | PLS_INTG_MDSplt1) /
    wrapped_legend +
    plot_layout(heights = c(4,4,2)) +
    plot_annotation(title = "MDS Plots: Outlining the Integrated Results for DNA and RNA with and without Time-Scale\n", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5),
                                  plot.margin = unit(c(0.5, 3, 0, 0), "cm")))



################################################################################
# Adding the PCA and PLS Bray-Curtis comparisons


PCA_DNA_MDSplt0
PCA_RNA_MDSplt0












