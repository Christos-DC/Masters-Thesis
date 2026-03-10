
# This is getting the MDS plots for the Bray-Curtis, Jaccard and Hill-based metrics
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")


# Defining the parameters
conditions <- DNA.3$condition
reactors <- DNA.3$reactor
timepts <- DNA.3$week_from_salt1


# Bray-Curtis (DNA)
DNABray <- sampledist(DNA.1, Bray_Curtis)
title <- "MDS plot for Bray-Curtis(DNA)"

MDSstructure <- mdsplotfunc(DNABray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))


numerics <- MDSnumeric(df = DNABray,
           func = PCA_Bray_Curtis,
           reactor = reactors,
           condition = conditions,
           timepts = timepts,
           distdf = TRUE)
numerics


# Bray-Curtis (RNA)
RNABray <- sampledist(RNA.1, Bray_Curtis)
title <- "MDS plot for Bray-Curtis (RNA)"

MDSstructure <- mdsplotfunc(RNABray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))



# Jaccard (DNA)
DNAJacc <- sampledist(DNA.1, Jaccard)
title <- "MDS plot for Jaccard (DNA)"

MDSstructure <- mdsplotfunc(DNAJacc, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))



# Jaccard (RNA)
RNAJacc <- sampledist(RNA.1, Jaccard)
title <- "MDS plot for Jaccard (RNA)"

MDSstructure <- mdsplotfunc(RNAJacc, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))



# Hill-based (DNA)
q <- 1
DNAHill <- sampledist(DNA.1, hill, q=q)
title <- paste("MDS plot for Hill-based (DNA) with q =", q)

MDSstructure <- mdsplotfunc(DNAHill, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))


# Hill-based (RNA)
q <- 1
RNAHill <- sampledist(RNA.1, hill, q=q)
title <- paste("MDS plot for Hill-based (RNA) with q =", q)

MDSstructure <- mdsplotfunc(RNAHill, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))






















