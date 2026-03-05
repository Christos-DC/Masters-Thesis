
# This is the N-integration Step for the Zymo Dataset
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")


################################################################################

# This is just a simple N-integration step without any time adjustment.
ncomps <- 4
X1 <- logratio.transfo(DNA.1 + 0.001, logratio = "CLR")
X2 <- logratio.transfo(RNA.1 + 0.001, logratio = "CLR")

pls.result <- pls(X=X1, Y=X2, ncomp = ncomps, mode = "canonical")

compsDNA <- pls.result$variates$X
compsRNA <- pls.result$variates$Y

# View(compsDNA)
# View(compsRNA)
rownames(compsRNA) <- rownames(RNA.1)


# Defining the parameters
conditions <- DNA.3$condition
reactors <- DNA.3$reactor
timepts <- DNA.3$week_from_salt1

# For the DNA component
DNA_PLSBray <- sampledist(compsDNA, PCA_Bray_Curtis)
title <- paste("MDS plot: PLS Bray-Curtis (DNA) with", ncomps, "components.", sep = " ")

MDSstructure <- mdsplotfunc(DNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))



# For the RNA component
RNA_PLSBray <- sampledist(compsRNA, PCA_Bray_Curtis)
title <- paste("MDS plot: PLS Bray-Curtis (RNA) with", ncomps, "components.", sep = " ")

MDSstructure <- mdsplotfunc(RNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))



# What does averaging the components do?
compsCombined <- (compsDNA + compsRNA)/2

Comb_PLSBray <- sampledist(compsCombined, PCA_Bray_Curtis)
title <- paste("MDS plot: PLS Bray-Curtis (DNA and RNA) with", ncomps, "components.", sep = " ")

MDSstructure <- mdsplotfunc(Comb_PLSBray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))


################################################################################

# Doing the same thing as above but now with a linear kernel and lambda = 1
ncomps <- 4
X1 <- logratio.transfo(DNA.1 + 0.001, logratio = "CLR")
X2 <- logratio.transfo(RNA.1 + 0.001, logratio = "CLR")

pls.result <- pls(X=X1, Y=X2, ncomp = ncomps, mode = "canonical")

compsDNA <- pls.result$variates$X
compsRNA <- pls.result$variates$Y

rownames(compsRNA) <- rownames(RNA.1)


# Defining the parameters
conditions <- DNA.3$condition
reactors <- DNA.3$reactor
timepts <- DNA.3$week_from_salt1

# For the DNA component
# DNA_PCABray <- sampledist(compsDNA, PCA_Bray_Curtis)
lambda <- 1
DNA_PLSBray <- penaltyfunc(df = compsDNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
title <- paste("MDS plot: PLS Bray-Curtis (DNA) with", ncomps, "components. Using Linear Kernel", sep = " ")

MDSstructure <- mdsplotfunc(DNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))



# For the RNA component
RNA_PLSBray <- penaltyfunc(df = compsRNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
title <- paste("MDS plot: PLS Bray-Curtis (RNA) with", ncomps, "components. Using Linear Kernel", sep = " ")

MDSstructure <- mdsplotfunc(RNA_PLSBray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))


# What does averaging the components do?
compsCombined <- (compsDNA + compsRNA)/2

Comb_PLSBray <- penaltyfunc(df = compsCombined,
                            metric = PCA_Bray_Curtis,
                            reactor = reactors,
                            timepts = timepts,
                            kernelfunc = linearkernel,
                            lambda = lambda)
title <- paste("MDS plot: PLS Bray-Curtis (DNA and RNA) with", ncomps, "components. Using Linear Kernel", sep = " ")

MDSstructure <- mdsplotfunc(Comb_PLSBray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))






################################################################################

# Comparing the PCA and PLS approaches (when Lambda = 0)
ncomps <- 4

# Defining the parameters
conditions <- DNA.3$condition
reactors <- DNA.3$reactor
timepts <- DNA.3$week_from_salt1

# Get PCA components
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
PCAcompsDNA <- pca.DNA$variates$X

pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
PCAcompsRNA <- pca.RNA$variates$X

# Get PLS components
X1 <- logratio.transfo(DNA.1 + 0.001, logratio = "CLR")
X2 <- logratio.transfo(RNA.1 + 0.001, logratio = "CLR")

pls.result <- pls(X=X1, Y=X2, ncomp = ncomps, mode = "canonical")

PLScompsDNA <- pls.result$variates$X
PLScompsRNA <- pls.result$variates$Y

rownames(compsRNA) <- rownames(RNA.1)


# Comparing DNA (lambda = 0)
DNA_PCABray <- sampledist(PCAcompsDNA, PCA_Bray_Curtis)
DNA_PLSBray <- sampledist(PLScompsDNA, PCA_Bray_Curtis)

MDS_DNA_PCA <- mdsplotfunc(DNA_PCABray, conditions, timepts, "PCA")
MDS_DNA_PLS <- mdsplotfunc(DNA_PLSBray, conditions, timepts, "PLS")

MDS_DNA_PCA$stress
MDS_DNA_PLS$stress

DNA_PCAplt <- MDS_DNA_PCA$mdsplot
DNA_PLSplt <- MDS_DNA_PLS$mdsplot

title <- paste("PCA vs PLS on DNA data using", ncomps, "components.", sep = " ")
combined_plt <- (DNA_PCAplt + theme(legend.position = "none") | DNA_PLSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5, 1.3)) +
    plot_annotation(title = title, 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))

combined_plt


# Comparing RNA (lambda = 0)
RNA_PCABray <- sampledist(PCAcompsRNA, PCA_Bray_Curtis)
RNA_PLSBray <- sampledist(PLScompsRNA, PCA_Bray_Curtis)

MDS_RNA_PCA <- mdsplotfunc(RNA_PCABray, conditions, timepts, "PCA")
MDS_RNA_PLS <- mdsplotfunc(RNA_PLSBray, conditions, timepts, "PLS")

MDS_RNA_PCA$stress
MDS_RNA_PLS$stress

RNA_PCAplt <- MDS_RNA_PCA$mdsplot
RNA_PLSplt <- MDS_RNA_PLS$mdsplot

title <- paste("PCA vs PLS on RNA data using", ncomps, "components.", sep = " ")
combined_plt <- (RNA_PCAplt + theme(legend.position = "none") | RNA_PLSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5, 1.3)) +
    plot_annotation(title = title, 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))

combined_plt




# Comparing DNA (lambda = 1)
lambda <- 1
DNA_PCABray <- penaltyfunc(df = PCAcompsDNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
DNA_PLSBray <- penaltyfunc(df = PLScompsDNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)

MDS_DNA_PCA <- mdsplotfunc(DNA_PCABray, conditions, timepts, "PCA")
MDS_DNA_PLS <- mdsplotfunc(DNA_PLSBray, conditions, timepts, "PLS")

MDS_DNA_PCA$stress
MDS_DNA_PLS$stress

DNA_PCAplt <- MDS_DNA_PCA$mdsplot
DNA_PLSplt <- MDS_DNA_PLS$mdsplot

title <- paste("PCA vs PLS on DNA data using", ncomps, "components.", sep = " ")
combined_plt <- (DNA_PCAplt + theme(legend.position = "none") | DNA_PLSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5, 1.3)) +
    plot_annotation(title = title, 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))

combined_plt



# Comparing RNA (lambda = 1)
lambda <- 1
RNA_PCABray <- penaltyfunc(df = PCAcompsRNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)
RNA_PLSBray <- penaltyfunc(df = PLScompsRNA,
                           metric = PCA_Bray_Curtis,
                           reactor = reactors,
                           timepts = timepts,
                           kernelfunc = linearkernel,
                           lambda = lambda)

MDS_RNA_PCA <- mdsplotfunc(RNA_PCABray, conditions, timepts, "PCA")
MDS_RNA_PLS <- mdsplotfunc(RNA_PLSBray, conditions, timepts, "PLS")

MDS_RNA_PCA$stress
MDS_RNA_PLS$stress

RNA_PCAplt <- MDS_RNA_PCA$mdsplot
RNA_PLSplt <- MDS_RNA_PLS$mdsplot

title <- paste("PCA vs PLS on RNA data using", ncomps, "components.", sep = " ")
combined_plt <- (RNA_PCAplt + theme(legend.position = "none") | RNA_PLSplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5, 1.3)) +
    plot_annotation(title = title, 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))

combined_plt














 









