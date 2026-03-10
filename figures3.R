
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
inhib1_time <- df[87:98, ]
inhib2_time <- df[12:20, ]

# For the Inhibition 1 group
# inhib1_time[13,] <- c("inhib_1", 230, 0)
# inhib1_time[,2] <- as.integer(inhib1_time[,2])
# inhib1_time[,3] <- as.numeric(inhib1_time[,3])

timepts <- seq(0, 446)
timelength <- length(timepts)
df1 <- matrix(0, nrow = timelength, ncol = 2)
colnames(df1) <- c("Day", "Salt")
df1[,1] <- timepts

n <- dim(inhib1_time)[1]
for (i in 2:n){
    pt1 <- inhib1_time[,2][i-1]
    pt2 <- inhib1_time[,2][i]
    salt1 <- inhib1_time[,3][i-1]
    salt2 <- inhib1_time[,3][i]
    
    bw_time <- seq(pt1, pt2)
    salt_pred <- rev(seq(salt2, salt1, length.out = pt2 - pt1 + 1))
    
    df1[bw_time + 1, 2] <- salt_pred
}

df1 <- as.data.frame(df1)
plt1 <- df1 %>% ggplot(aes(x=Day, y = Salt)) +
    geom_point() +
    theme_minimal() +
    annotate("text", x = 100, y = 25, label = "Inhibition 1", colour = "blue") +
    geom_rect(aes(xmin = 125, xmax = 200, ymin = -0.1, ymax = 30), fill = NA, color = "blue",
              linewidth = 1) 


# For the Control and Inhibition 2 group
df2 <- matrix(0, nrow = timelength, ncol = 2)
colnames(df2) <- c("Day", "Salt")
df2[,1] <- timepts

n <- dim(inhib2_time)[1]
for (i in 2:n){
    pt1 <- inhib2_time[,2][i-1]
    pt2 <- inhib2_time[,2][i]
    salt1 <- inhib2_time[,3][i-1]
    salt2 <- inhib2_time[,3][i]
    
    bw_time <- seq(pt1, pt2)
    print(bw_time)
    salt_pred <- rev(seq(salt2, salt1, length.out = pt2 - pt1 + 1))
    
    df2[bw_time + 1, 2] <- salt_pred
}

df2 <- as.data.frame(df2)
plt2 <- df2 %>% ggplot(aes(x=Day, y = Salt)) +
    geom_point() +
    theme_minimal() +
    annotate("text", x = 105, y = 25, label = "Control", colour = "red") +
    geom_rect(aes(xmin = 125, xmax = 200, ymin = -0.1, ymax = 30), fill = NA, color = "red",
              linewidth = 1) +
    annotate("text", x = 355, y = 25, label = "Inhibition 2", colour = "green") +
    geom_rect(aes(xmin = 380, xmax = 440, ymin = -0.1, ymax = 30), fill = NA, color = "green",
              linewidth = 1) 
    

(plt2 / plt1) +
    plot_layout(axis_titles = "collect") +
    plot_annotation(title = "Salt concentration in the CSTR", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)))


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




































