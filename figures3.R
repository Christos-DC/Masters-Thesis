
#' This is a file that would be used to outline the figures used for my thesis
#' that are made in R.
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")

################################################################################

# Figure 3.2 - Salt Concentration
which(DNA.3$NaCl..g.L. != RNA.3$NaCl..g.L.)

SaltCon <- DNA.3$NaCl..g.L. 
day <- DNA.3$day
df <- data.frame("condition" = conditions, "Day" = day, "Salt" = SaltCon)

# use the sequence function to get values in between perhaps?
control_time <- df[1:9, ]
inhib1_time <- df[87:98, ]
inhib2_time <- df[12:20, ]

# For the Inhibition 1 group
timepts <- seq(0, 446)
timelength <- length(timepts)
df1 <- matrix(0, nrow = timelength, ncol = 2)
colnames(df1) <- c("Day", "Salt")
df1[,1] <- timepts

# Outlining the general trend for salt concentration in experiment
expinterp <- function(x, pt1, pt2){
    x1 <- pt1[1]
    x2 <- pt2[1]
    y1 <- pt1[2]
    y2 <- pt2[2]
    
    b <- (y2/y1)^(1/(x2-x1))
    a <- y1 / b^x1
    
    return(a*b^x)
}

pt1 <- as.numeric(inhib1_time[1, c(2,3)])
pt2 <- as.numeric(inhib1_time[12, c(2,3)])
df1[131:299 + 1, 2] <- expinterp(131:299, pt1, pt2)


df1 <- as.data.frame(df1)
plt1 <- df1 %>% ggplot(aes(x=Day, y = Salt)) +
    geom_point(size = 2) +
    theme_minimal() +
    annotate("text", x = 75, y = 25, label = "Inhibition 1", colour = "blue", size = 6) +
    geom_rect(aes(xmin = 125, xmax = 210, ymin = -0.1, ymax = 30), fill = NA, color = "blue",
              linewidth = 1) +
    theme(axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "Salt Concentration")


# For the Control and Inhibition 2 group
df2 <- matrix(0, nrow = timelength, ncol = 2)
colnames(df2) <- c("Day", "Salt")
df2[,1] <- timepts


pt1 <- as.numeric(inhib2_time[1, c(2,3)])
pt2 <- as.numeric(inhib2_time[9, c(2,3)])
df2[385:446 + 1, 2] <- expinterp(385:446, pt1, pt2)


df2 <- as.data.frame(df2)
plt2 <- df2 %>% ggplot(aes(x=Day, y = Salt)) +
    geom_point(size = 2) +
    theme_minimal() +
    annotate("text", x = 90, y = 25, label = "Control", colour = "red", size = 6) +
    geom_rect(aes(xmin = 125, xmax = 200, ymin = -0.1, ymax = 30), fill = NA, color = "red",
              linewidth = 1) +
    annotate("text", x = 330, y = 25, label = "Inhibition 2", colour = "green", size = 6) +
    geom_rect(aes(xmin = 380, xmax = 450, ymin = -0.1, ymax = 30), fill = NA, color = "green",
              linewidth = 1) +
    theme(axis.title = element_text(size = 16),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "Salt Concentration")
    
# Combined Plot
(plt2 / plt1) +
    plot_layout(axis_titles = "collect") +
    plot_annotation(title = "Salt concentration in the Continuous Stirred-Tank Reactor", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))


################################################################################

# Figure 3.3 - Bray-Curtis on DNA and RNA data

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

# Figure 3.4 - Jaccard on DNA and RNA data.

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


MDS_DNAplt <- MDS_DNAplt + theme(legend.position = "none")
MDS_RNAplt <- MDS_RNAplt + theme(axis.title.y = element_blank())

plt <- (MDS_DNAplt | MDS_RNAplt) / wrapped_legend + 
    plot_layout(heights = c(5,1.5)) +
    plot_annotation(title = "NMDS on Jaccard Distance Matrix", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))
plt

ggsave("Jaccard.png", plot = plt, path = "./figures/",
       width = 24, height = 14, units = "cm")



################################################################################

# Figure - Hill-based on DNA and RNA data for q=0,1,2

title <- "DNA"
DNAHill0 <- sampledist(DNA.1, hill, q = 0)
DNAHill1 <- sampledist(DNA.1, hill, q = 1)
DNAHill2 <- sampledist(DNA.1, hill, q = 2)

MDSstructure <- mdsplotfunc(DNAHill0, conditions, timepts, title)
MDSstructure$stress
MDS_DNAplt0 <- MDSstructure$mdsplot

MDSstructure <- mdsplotfunc(DNAHill1, conditions, timepts, title)
MDSstructure$stress
MDS_DNAplt1 <- MDSstructure$mdsplot

MDSstructure <- mdsplotfunc(DNAHill2, conditions, timepts, title)
MDSstructure$stress
MDS_DNAplt2 <- MDSstructure$mdsplot


title <- "RNA"
RNAHill0 <- sampledist(RNA.1, hill, q = 0)
RNAHill1 <- sampledist(RNA.1, hill, q = 1)
RNAHill2 <- sampledist(RNA.1, hill, q = 2)

MDSstructure <- mdsplotfunc(RNAHill0, conditions, timepts, title)
MDSstructure$stress
MDS_RNAplt0 <- MDSstructure$mdsplot

MDSstructure <- mdsplotfunc(RNAHill1, conditions, timepts, title)
MDSstructure$stress
MDS_RNAplt1 <- MDSstructure$mdsplot

MDSstructure <- mdsplotfunc(RNAHill2, conditions, timepts, title)
MDSstructure$stress
MDS_RNAplt2 <- MDSstructure$mdsplot


# Adjusting the plots
MDS_DNAplt0 <- MDS_DNAplt0 + theme(axis.title.x = element_blank(), 
                                   legend.position = "none") +
    ylab(expression(atop("q=0", "D2"))) 

MDS_DNAplt1 <- MDS_DNAplt1 + theme(axis.title.x = element_blank(), 
                                   legend.position = "none",
                                   plot.title = element_blank()) +
    ylab(expression(atop("q=1", "D2")))

MDS_DNAplt2 <- MDS_DNAplt2 + theme(legend.position = "none",
                                   plot.title = element_blank()) +
    ylab(expression(atop("q=2", "D2")))


MDS_RNAplt0 <- MDS_RNAplt0 + theme(axis.title.x = element_blank(), 
                                   axis.title.y = element_blank(),
                                   legend.position = "none") 

MDS_RNAplt1 <- MDS_RNAplt1 + theme(axis.title.x = element_blank(), 
                                   axis.title.y = element_blank(),
                                   plot.title = element_blank()) 

MDS_RNAplt2 <- MDS_RNAplt2 + theme(axis.title.y = element_blank(),
                                   legend.position = "none",
                                   plot.title = element_blank()) 


# Combined plot
plt <- (MDS_DNAplt0 | MDS_RNAplt0) /
    (MDS_DNAplt1 | MDS_RNAplt1) /
    (MDS_DNAplt2 | MDS_RNAplt2) / 
    wrapped_legend +
    plot_layout(heights = c(4,4,4,1.5)) +
    plot_annotation(title = "NMDS on Hill-based Distance Matrix for q=0,1,2", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))

plt

ggsave("Hill.png", plot = plt, path = "./figures/",
       width = 24, height = 28, units = "cm")



################################################################################























