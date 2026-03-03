
# This is the N-integration Step for the Zymo Dataset
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")


########### Learning how to perform N-integration on a simple dataset.##########

# This is how the design matrix is set up
design = matrix(1, ncol = 3, nrow = 3, 
                dimnames = list(c("mirna", "mrna", "protein"), 
                                c("mirna", "mrna", "protein")))
diag(design) = 0
design


# This is what is happening with actual datasets.
data(breast.TCGA) # extract the TCGA data

X1 <- breast.TCGA$data.train$mirna # use the mirna and mrna expression levels as 
X2 <- breast.TCGA$data.train$mrna  # the X datasets
X <- list(mirna = X1, mrna = X2)

Y <- breast.TCGA$data.train$protein # set the protein levels as the Y dataset

block.pls.result <- block.pls(X, Y, design = "full") # run the method

plotIndiv(block.pls.result) # plot the samples
plotVar(block.pls.result, legend = TRUE) # plot the variables



################################################################################
ncomps <- 4
pls.result <- pls(RNA.1 + 0.001, DNA.1 + 0.001, ncomp = ncomps, mode = "canonical", logratio = 'CLR')

# View(pls.result$variates$X)
# View(pls.result$variates$Y)

compsDNA <- pls.result$variates$Y
compsRNA <- pls.result$variates$X

# Defining the parameters
conditions <- DNA.3$condition
reactors <- DNA.3$reactor
timepts <- DNA.3$week_from_salt1

# For the DNA component
DNA_PCABray <- sampledist(compsDNA, PCA_Bray_Curtis)
title <- paste("MDS plot: PLS Bray-Curtis (DNA) with", ncomps, "components.", sep = " ")

MDSstructure <- mdsplotfunc(DNA_PCABray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))



# For the RNA component
RNA_PCABray <- sampledist(compsRNA, PCA_Bray_Curtis)
title <- paste("MDS plot: PLS Bray-Curtis (RNA) with", ncomps, "components.", sep = " ")

MDSstructure <- mdsplotfunc(RNA_PCABray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))


#' The problem is symmetrical. Whichever dataset is set in the Y component, we 
#' don't have a nice plot that is similar to what is shown in PCA. However the 
#' X dataset does behave in a similar manner. I am not sure why but this needs 
#' to be addressed.
#' Is it due to PLS having to use Y as the response variable and it's doing something
#' else underneath the function call?




 









