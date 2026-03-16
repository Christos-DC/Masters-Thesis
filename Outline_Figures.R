

# Figures for Oliver to see
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")

# Defining the parameters
conditions <- DNA.3$condition
reactors <- DNA.3$reactor
timepts <- DNA.3$week_from_salt1

################################################################################
# My PCA plots
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = "CLR")
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = "CLR")

# Plotting to check the correlation between the two datasets.
hex_colors <- blend_with_white(base_colours[conditions], timepts)

compsDNA <- pca.DNA$variates$X
compsRNA <- pca.RNA$variates$X


# Making my own PLS plots for DNA and RNA
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
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(size = 12, hjust = 0.5))

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
          plot.title = element_text(size = 12, hjust = 0.5))

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

# This is just a simple N-integration step without any time adjustment.
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
          panel.background = element_rect(fill = "white", color = NA),
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
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(size = 12, hjust = 0.5))

(pls_plt1 + theme(legend.position = "none") | 
        pls_plt2 + theme(axis.title.y = element_blank())) / wrapped_legend +
    plot_layout(heights = c(5,1.3)) +
    plot_annotation(title = "PLS plots for DNA and RNA", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))






# No Time adjustment



























