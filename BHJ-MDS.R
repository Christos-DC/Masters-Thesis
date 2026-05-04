
# This is getting the MDS plots for the Bray-Curtis, Jaccard and Hill-based metrics
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")


# Bray-Curtis (DNA)
DNABray <- sampledist(DNA.1, BrayCurtis)
title <- "MDS plot for Bray-Curtis(DNA)"

MDSstructure <- mdsplotfunc(DNABray, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))


numerics <- MDSnumeric(df = DNABray,
           func = BrayCurtis,
           reactor = reactors,
           condition = conditions,
           timepts = timepts,
           distdf = TRUE)
numerics


# Bray-Curtis (RNA)
RNABray <- sampledist(RNA.1, BrayCurtis)
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
q <- 2
DNAHill <- sampledist(DNA.1, Hill, q=q)
title <- paste("MDS plot for Hill-based (DNA) with q =", q)

MDSstructure <- mdsplotfunc(DNAHill, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))


# Hill-based (RNA)
q <- 2
RNAHill <- sampledist(RNA.1, Hill, q=q)
title <- paste("MDS plot for Hill-based (RNA) with q =", q)

MDSstructure <- mdsplotfunc(RNAHill, conditions, timepts, title)
MDSstructure$stress

MDSplt <- MDSstructure$mdsplot
MDSplt / wrapped_legend + plot_layout(heights = c(5, 1.3))












