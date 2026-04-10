
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


































