
#' This is a file that would be used to outline the figures used for my thesis
#' that are made in R.
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")

# Defining the parameters
conditions <- DNA.3$condition
reactors <- DNA.3$reactor
timepts <- DNA.3$week_from_salt1

################################################################################

# Replicating the Salt concentration figure from Olivier
which(DNA.3$NaCl..g.L. != RNA.3$NaCl..g.L.)

SaltCon <- DNA.3$NaCl..g.L. 
day <- DNA.3$day
df <- data.frame("condition" = conditions, "Day" = day, "Salt" = SaltCon)

# use the sequence function to get values in between perhaps?
control_time <- df[1:9, ]
inhib1_time <- df[85:98, ]
inhib2_time <- df[10:20, ]

timepts <- seq(0, 450)
inhib1_time





################################################################################

# Figure 3.2 - Bray-Curtis on DNA and RNA data

DNABray <- sampledist(DNA.1, Bray_Curtis)
title <- "DNA"

MDSstructure <- mdsplotfunc(DNABray, conditions, timepts, title)
MDSstructure$stress

MDS_DNAplt <- MDSstructure$mdsplot


RNABray <- sampledist(RNA.1, Bray_Curtis)
title <- "RNA"

MDSstructure <- mdsplotfunc(RNABray, conditions, timepts, title)
MDSstructure$stress

MDS_RNAplt <- MDSstructure$mdsplot

plt <- (MDS_DNAplt + theme(legend.position = "none")| MDS_RNAplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
            plot_layout(heights = c(5,1.5)) +
    plot_annotation(title = "NMDS on Bray-Curtis Distance Matrix", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))

ggsave("Bray-Curtis.png", plot = plt, path = "./figures/",
       width = 24, height = 14, units = "cm")

# ggsave("Bray-Curtis.png", plot = plt, path = "./figures/",
#        width = 18, height = 20, units = "cm")

################################################################################

# Figure 3.3 - Jaccard on DNA and RNA data.

DNAJacc <- sampledist(DNA.1, Jaccard)
title <- "DNA"

MDSstructure <- mdsplotfunc(DNAJacc, conditions, timepts, title)
MDSstructure$stress

MDS_DNAplt <- MDSstructure$mdsplot


RNAJacc <- sampledist(RNA.1, Jaccard)
title <- "RNA"

MDSstructure <- mdsplotfunc(RNAJacc, conditions, timepts, title)
MDSstructure$stress

MDS_RNAplt <- MDSstructure$mdsplot

plt <- (MDS_DNAplt + theme(legend.position = "none")| MDS_RNAplt + theme(axis.title.y = element_blank())) / wrapped_legend + 
    plot_layout(heights = c(5,1.5)) +
    plot_annotation(title = "NMDS on Jaccard Distance Matrix", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))
plt

ggsave("Jaccard.png", plot = plt, path = "./figures/",
       width = 24, height = 14, units = "cm")




































