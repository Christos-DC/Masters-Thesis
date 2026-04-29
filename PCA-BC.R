
# This is testing for the PCA Bray-Curtis metric and the Time-adjustment metric
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")


# PCA Bray-Curtis (DNA)
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
comps <- pca.DNA$variates$X

# Find out the meanings of the components used for plotting
# plotLoadings(pca.DNA, comp = 1, contrib = 'max', method = 'mean')
# plotLoadings(pca.DNA, comp = 2, contrib = 'max', method = 'mean')

# top_vars_comp1 <- selectVar(pca.DNA, comp = 1)$value
# head(top_vars_comp1)

DNA_PCABray <- sampledist(comps, PCA_Bray_Curtis)
title <- paste("MDS plot: PCA Bray-Curtis (DNA) with", ncomps, "components.", sep = " ")

MDSstructure <- mdsplotfunc(DNA_PCABray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))


# PCA Bray-Curtis (RNA)
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.RNA$variates$X[, 1:ncomps]

RNA_PCABray <- sampledist(comps, PCA_Bray_Curtis)
title <- paste("MDS plot: PCA Bray-Curtis (RNA) with", ncomps, "components.", sep = " ")

MDSstructure <- mdsplotfunc(RNA_PCABray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))



# Time-adjusted Bray-Curtis (DNA)
lambdas <- seq(0, 1, by = 0.05)
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.DNA$variates$X[, 1:ncomps]


dfs <- list() #this holds all the distance matrices that have been adjusted based on a particular lambda value.
k <- 1
for (lambda in lambdas){
    dfs[[k]] <- penaltyfunc(comps, PCA_Bray_Curtis, reactors, timepts, linearkernel, lambda)
    k <- k + 1
}

mdsplots <- list()
for (i in 1:length(lambdas)){
    df <- dfs[[i]]
    df <- df[,-c(99,100)]
    lambda <- lambdas[i]
    text <- paste("MDS plot: PCA Bray-Curtis (DNA) with ", ncomps, " components, Lambda = ", lambda,", Linear Kernel", sep = "")
    mdsplots[[i]] <- mdsplotfunc(df,conditions, timepts, text)
}


mdsplots[[21]]



# Time-adjusted Bray-Curtis (RNA)
lambdas <- seq(0, 1, by = 0.05)
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.RNA$variates$X[, 1:ncomps]


dfs <- list() #this holds all the distance matrices that have been adjusted based on a particular lambda value.
k <- 1
for (lambda in lambdas){
    dfs[[k]] <- penaltyfunc(comps, PCA_Bray_Curtis, reactors, timepts, linearkernel, lambda)
    k <- k + 1
}

mdsplots <- list()
for (i in 1:length(lambdas)){
    df <- dfs[[i]]
    df <- df[,-c(99,100)]
    lambda <- lambdas[i]
    text <- paste("MDS plot: PCA Bray-Curtis (RNA) with ", ncomps, " components, Lambda = ", lambda,", Linear Kernel", sep = "")
    mdsplots[[i]] <- mdsplotfunc(df,conditions, timepts, text)
}


mdsplots[[1]]



################################################################################

# Analysing the stress value over time
n <- length(lambdas)
stressvals <- rep(0, n)
for (i in 1:n){
    stressvals[i] <- mdsplots[[i]]$stress
}

plot(lambdas, stressvals, type = 'l', xlab = "λ values", ylab = "Stress Values",
     main = "Plotting stress values from the MDS plots for corresponding λ")


#' The stress values tend to be low enough for it to be an acceptable plot
#' but it's high enough where we have to have some cautions on it.
#' This is likely because I am using NMDS which is not using a true metric.
#' Possible way to get around this is by either using or developing a new 
#' true metric that can possibly make the mds results have a lower stress value.

# Since I am performing mds on the PCA data transformation, why don't I try the euclidean metric.

ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = 10, logratio = 'CLR')
comps <- pca.DNA$variates$X[, 1:ncomps]

DNA_PCABray <- sampledist(comps, euclidean)
title <- paste("MDS plot: PCA Bray-Curtis (DNA) with", ncomps, "components.", sep = " ")

MDSstructure <- mdsplotfunc(DNA_PCABray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))

#' Plots are not the most interpretable as expected but it supports the idea that
#' the stress value is now lower since I am using a valid metric instead of a 
#' semi-metric. 











