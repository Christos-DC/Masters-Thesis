
# This is the code for the second appendix chapter I have
source("~/Documents/Masters Degree/Masters Research/Code Scripts/Polished Code/data.R")


##################### DNA - Linear kernel function #############################
lambdas <- seq(0, 1, by = 0.05)
n <- length(lambdas)
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
comps <- pca.DNA$variates$X

# Setting up the matrices needed for plotting (DNA)
lambdadf <- matrix(NA, nrow = n, ncol = 4)
colnames(lambdadf) <- c("lambda", "control", "inhib_1", "inhib_2")
CI_vals <- matrix(NA, nrow = n, ncol = 7)
colnames(CI_vals) <- c("lambda", "C_low", "C_high", "I1_low", "I1_high", "I2_low", "I2_high")


# Gathering the numerical information 
for (i in 1:n){
    lambda <- lambdas[i]
    df <- penaltyfunc(comps, CompBC, reactors, timepts, linearkernel, lambda)
    DNAnumerics <- MDSnumeric(df = df[, -c(99,100)], 
                              func = CompBC, 
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
    
    meanvals <- DNAnumerics$meanvals
    CI <- DNAnumerics$confint
    lambdadf[i, ] <- c(lambda, meanvals)
    CI_vals[i,] <- c(lambda, CI[1,], CI[2,], CI[3,])
    
    if (i %% 10 == 0) print(i)
}

# Plotting the ribbon plots
lambdadf <- as.data.frame(lambdadf)
base_colors <- c("control" = "red", "inhib_1" = "blue", "inhib_2" = "green")

LinearDNA <- lambdadf %>% ggplot(aes(x=lambda)) +
    geom_line(aes(y = control, col = "control"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = C_low, ymax= C_high), fill = "red", alpha = 0.2) +
    geom_line(aes(y = inhib_1, col = "inhib_1"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I1_low, ymax= I1_high), fill = "blue", alpha = 0.2) +
    geom_line(aes(y = inhib_2, col = "inhib_2"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I2_low, ymax= I2_high), fill = "green", alpha = 0.2) +
    theme_minimal() +
    labs(title = expression(bold("DNA")), 
         y = expression(atop(bold("Linear"), "Mean Distance Values")),
         col = "") +
    scale_color_manual(values = base_colors) +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 12, hjust = 0.5))



##################### RNA - Linear kernel function #############################
lambdas <- seq(0, 1, by = 0.05)
n <- length(lambdas)
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
comps <- pca.RNA$variates$X

# Setting up the matrices needed for plotting (DNA)
lambdadf <- matrix(NA, nrow = n, ncol = 4)
colnames(lambdadf) <- c("lambda", "control", "inhib_1", "inhib_2")
CI_vals <- matrix(NA, nrow = n, ncol = 7)
colnames(CI_vals) <- c("lambda", "C_low", "C_high", "I1_low", "I1_high", "I2_low", "I2_high")


# Gathering the numerical information 
for (i in 1:n){
    lambda <- lambdas[i]
    df <- penaltyfunc(comps, CompBC, reactors, timepts, linearkernel, lambda)
    DNAnumerics <- MDSnumeric(df = df[, -c(99,100)], 
                              func = CompBC, 
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
    
    meanvals <- DNAnumerics$meanvals
    CI <- DNAnumerics$confint
    lambdadf[i, ] <- c(lambda, meanvals)
    CI_vals[i,] <- c(lambda, CI[1,], CI[2,], CI[3,])
    
    if (i %% 10 == 0) print(i)
}

# Plotting the ribbon plots
lambdadf <- as.data.frame(lambdadf)
base_colors <- c("control" = "red", "inhib_1" = "blue", "inhib_2" = "green")

LinearRNA <- lambdadf %>% ggplot(aes(x=lambda)) +
    geom_line(aes(y = control, col = "control"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = C_low, ymax= C_high), fill = "red", alpha = 0.2) +
    geom_line(aes(y = inhib_1, col = "inhib_1"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I1_low, ymax= I1_high), fill = "blue", alpha = 0.2) +
    geom_line(aes(y = inhib_2, col = "inhib_2"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I2_low, ymax= I2_high), fill = "green", alpha = 0.2) +
    theme_minimal() +
    labs(title = expression(bold("RNA")), 
         col = "") +
    scale_color_manual(values = base_colors) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 12, hjust = 0.5))



##################### DNA - Exponential kernel function ########################
lambdas <- seq(0, 1, by = 0.05)
n <- length(lambdas)
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
comps <- pca.DNA$variates$X

# Setting up the matrices needed for plotting (DNA)
lambdadf <- matrix(NA, nrow = n, ncol = 4)
colnames(lambdadf) <- c("lambda", "control", "inhib_1", "inhib_2")
CI_vals <- matrix(NA, nrow = n, ncol = 7)
colnames(CI_vals) <- c("lambda", "C_low", "C_high", "I1_low", "I1_high", "I2_low", "I2_high")


# Gathering the numerical information 
for (i in 1:n){
    lambda <- lambdas[i]
    df <- penaltyfunc(comps, CompBC, reactors, timepts, expkernel, lambda)
    DNAnumerics <- MDSnumeric(df = df[, -c(99,100)], 
                              func = CompBC, 
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
    
    meanvals <- DNAnumerics$meanvals
    CI <- DNAnumerics$confint
    lambdadf[i, ] <- c(lambda, meanvals)
    CI_vals[i,] <- c(lambda, CI[1,], CI[2,], CI[3,])
    
    if (i %% 10 == 0) print(i)
}

# Plotting the ribbon plots
lambdadf <- as.data.frame(lambdadf)
base_colors <- c("control" = "red", "inhib_1" = "blue", "inhib_2" = "green")

ExpDNA <- lambdadf %>% ggplot(aes(x=lambda)) +
    geom_line(aes(y = control, col = "control"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = C_low, ymax= C_high), fill = "red", alpha = 0.2) +
    geom_line(aes(y = inhib_1, col = "inhib_1"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I1_low, ymax= I1_high), fill = "blue", alpha = 0.2) +
    geom_line(aes(y = inhib_2, col = "inhib_2"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I2_low, ymax= I2_high), fill = "green", alpha = 0.2) +
    theme_minimal() +
    labs(y = expression(atop(bold("Exponential"), "Mean Distance Values")),
         col = "") +
    scale_color_manual(values = base_colors) +
    theme(axis.title.x = element_blank(),
          legend.position = "none")


##################### RNA - Exponential kernel function ########################
lambdas <- seq(0, 1, by = 0.05)
n <- length(lambdas)
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
comps <- pca.RNA$variates$X

# Setting up the matrices needed for plotting (DNA)
lambdadf <- matrix(NA, nrow = n, ncol = 4)
colnames(lambdadf) <- c("lambda", "control", "inhib_1", "inhib_2")
CI_vals <- matrix(NA, nrow = n, ncol = 7)
colnames(CI_vals) <- c("lambda", "C_low", "C_high", "I1_low", "I1_high", "I2_low", "I2_high")


# Gathering the numerical information 
for (i in 1:n){
    lambda <- lambdas[i]
    df <- penaltyfunc(comps, CompBC, reactors, timepts, expkernel, lambda)
    DNAnumerics <- MDSnumeric(df = df[, -c(99,100)], 
                              func = CompBC, 
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
    
    meanvals <- DNAnumerics$meanvals
    CI <- DNAnumerics$confint
    lambdadf[i, ] <- c(lambda, meanvals)
    CI_vals[i,] <- c(lambda, CI[1,], CI[2,], CI[3,])
    
    if (i %% 10 == 0) print(i)
}

# Plotting the ribbon plots
lambdadf <- as.data.frame(lambdadf)
base_colors <- c("control" = "red", "inhib_1" = "blue", "inhib_2" = "green")

ExpRNA <- lambdadf %>% ggplot(aes(x=lambda)) +
    geom_line(aes(y = control, col = "control"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = C_low, ymax= C_high), fill = "red", alpha = 0.2) +
    geom_line(aes(y = inhib_1, col = "inhib_1"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I1_low, ymax= I1_high), fill = "blue", alpha = 0.2) +
    geom_line(aes(y = inhib_2, col = "inhib_2"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I2_low, ymax= I2_high), fill = "green", alpha = 0.2) +
    theme_minimal() +
    labs(col = "Condition") +
    scale_color_manual(values = base_colors) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text())




##################### DNA - Logarithmic kernel function ########################
lambdas <- seq(0, 1, by = 0.05)
n <- length(lambdas)
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
comps <- pca.DNA$variates$X

# Setting up the matrices needed for plotting (DNA)
lambdadf <- matrix(NA, nrow = n, ncol = 4)
colnames(lambdadf) <- c("lambda", "control", "inhib_1", "inhib_2")
CI_vals <- matrix(NA, nrow = n, ncol = 7)
colnames(CI_vals) <- c("lambda", "C_low", "C_high", "I1_low", "I1_high", "I2_low", "I2_high")


# Gathering the numerical information 
for (i in 1:n){
    lambda <- lambdas[i]
    df <- penaltyfunc(comps, CompBC, reactors, timepts, logkernel, lambda)
    DNAnumerics <- MDSnumeric(df = df[, -c(99,100)], 
                              func = CompBC, 
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
    
    meanvals <- DNAnumerics$meanvals
    CI <- DNAnumerics$confint
    lambdadf[i, ] <- c(lambda, meanvals)
    CI_vals[i,] <- c(lambda, CI[1,], CI[2,], CI[3,])
    
    if (i %% 10 == 0) print(i)
}

# Plotting the ribbon plots
lambdadf <- as.data.frame(lambdadf)
base_colors <- c("control" = "red", "inhib_1" = "blue", "inhib_2" = "green")

LogDNA <- lambdadf %>% ggplot(aes(x=lambda)) +
    geom_line(aes(y = control, col = "control"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = C_low, ymax= C_high), fill = "red", alpha = 0.2) +
    geom_line(aes(y = inhib_1, col = "inhib_1"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I1_low, ymax= I1_high), fill = "blue", alpha = 0.2) +
    geom_line(aes(y = inhib_2, col = "inhib_2"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I2_low, ymax= I2_high), fill = "green", alpha = 0.2) +
    theme_minimal() +
    labs(x = "λ Values",
         y = expression(atop(bold("Logarithmic"), "Mean Distance Values")),
         col = "") +
    scale_color_manual(values = base_colors) +
    theme(legend.position = "none")


##################### RNA - Exponential kernel function ########################
lambdas <- seq(0, 1, by = 0.05)
n <- length(lambdas)
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
comps <- pca.RNA$variates$X

# Setting up the matrices needed for plotting (DNA)
lambdadf <- matrix(NA, nrow = n, ncol = 4)
colnames(lambdadf) <- c("lambda", "control", "inhib_1", "inhib_2")
CI_vals <- matrix(NA, nrow = n, ncol = 7)
colnames(CI_vals) <- c("lambda", "C_low", "C_high", "I1_low", "I1_high", "I2_low", "I2_high")


# Gathering the numerical information 
for (i in 1:n){
    lambda <- lambdas[i]
    df <- penaltyfunc(comps, CompBC, reactors, timepts, logkernel, lambda)
    DNAnumerics <- MDSnumeric(df = df[, -c(99,100)], 
                              func = CompBC, 
                              reactor = reactors,
                              condition = conditions,
                              timepts = timepts,
                              distdf = TRUE)
    
    meanvals <- DNAnumerics$meanvals
    CI <- DNAnumerics$confint
    lambdadf[i, ] <- c(lambda, meanvals)
    CI_vals[i,] <- c(lambda, CI[1,], CI[2,], CI[3,])
    
    if (i %% 10 == 0) print(i)
}

# Plotting the ribbon plots
lambdadf <- as.data.frame(lambdadf)
base_colors <- c("control" = "red", "inhib_1" = "blue", "inhib_2" = "green")

LogRNA <- lambdadf %>% ggplot(aes(x=lambda)) +
    geom_line(aes(y = control, col = "control"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = C_low, ymax= C_high), fill = "red", alpha = 0.2) +
    geom_line(aes(y = inhib_1, col = "inhib_1"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I1_low, ymax= I1_high), fill = "blue", alpha = 0.2) +
    geom_line(aes(y = inhib_2, col = "inhib_2"), linewidth = 1) +
    geom_ribbon(data = CI_vals, aes(x=lambda, ymin = I2_low, ymax= I2_high), fill = "green", alpha = 0.2) +
    theme_minimal() +
    labs(x = "λ Values",
         col = "") +
    scale_color_manual(values = base_colors) +
    theme(axis.title.y = element_blank(),
          legend.position = "none")



############################# Plot the λ plots #################################

plt <- (LinearDNA | LinearRNA) / 
    (ExpDNA | ExpRNA) /
    (LogDNA | LogRNA) + 
    plot_layout(heights = c(4,4,4))

ggsave("lambdaChange.png", plot = plt, path = "./figures/",
       width = 24, height = 26, units = "cm")


################################################################################

# Here using DNA to plot the NMDS plots for λ=0,1 for lienar, exponential and log
ncomps <- 4
pca.DNA <- pca(DNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
DNAcomps <- pca.DNA$variates$X


# λ = 0
lambda <- 0
NoTime <- penaltyfunc(df = DNAcomps,
                            metric = CompBC,
                            reactor = reactors,
                            timepts = timepts,
                            kernelfunc = linearkernel,
                            lambda = lambda)
title <- expression(bold("No Reactor Time-Scale"))

MDSstructure <- mdsplotfunc(NoTime , conditions, timepts, title)
stress <- paste(100*round(MDSstructure$stress, 4), "%", sep = "")

MDS_NoTime  <- MDSstructure$mdsplot
MDS_NoTime  <- MDS_NoTime  + 
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    annotate("text", x = -Inf, y = -Inf, 
             label = stress, hjust = -0.15, vjust = -0.8,
             size = 5)


# Linear, λ = 1
lambda <- 1
Linear1 <- penaltyfunc(df = DNAcomps,
                       metric = CompBC,
                       reactor = reactors,
                       timepts = timepts,
                       kernelfunc = linearkernel,
                       lambda = lambda)
title <- "Linear Kernel"

MDSstructure <- mdsplotfunc(Linear1, conditions, timepts, title)
stress <- paste(100*round(MDSstructure$stress, 4), "%", sep = "")

MDS_Linear1 <- MDSstructure$mdsplot
MDS_Linear1 <- MDS_Linear1 + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    annotate("text", x = -Inf, y = -Inf, 
             label = stress, hjust = -0.15, vjust = -0.8,
             size = 5)



# Exponential, λ = 1
lambda <- 1
Exp1 <- penaltyfunc(df = DNAcomps,
                       metric = CompBC,
                       reactor = reactors,
                       timepts = timepts,
                       kernelfunc = expkernel,
                       lambda = lambda)
title <- "Exponential Kernel"

MDSstructure <- mdsplotfunc(Exp1, conditions, timepts, title)
stress <- paste(100*round(MDSstructure$stress, 4), "%", sep = "")

MDS_Exp1 <- MDSstructure$mdsplot
MDS_Exp1 <- MDS_Exp1 + 
    theme(legend.position = "none") +
    annotate("text", x = -Inf, y = -Inf, 
             label = stress, hjust = -0.15, vjust = -0.8,
             size = 5)


# Logarithmic, λ = 1
lambda <- 1
Log1 <- penaltyfunc(df = DNAcomps,
                    metric = CompBC,
                    reactor = reactors,
                    timepts = timepts,
                    kernelfunc = logkernel,
                    lambda = lambda)
title <- "Logarithmic Kernel"

MDSstructure <- mdsplotfunc(Log1, conditions, timepts, title)
stress <- paste(100*round(MDSstructure$stress, 4), "%", sep = "")

MDS_Log1 <- MDSstructure$mdsplot
MDS_Log1 <- MDS_Log1 + 
    theme(legend.position = c(1.15,1.1),
          axis.title.y = element_blank()) +
    annotate("text", x = -Inf, y = -Inf, 
             label = stress, hjust = -0.15, vjust = -0.8,
             size = 5)


plt <- (MDS_NoTime | MDS_Linear1) /
    (MDS_Exp1 | MDS_Log1) /
    wrapped_legend +
    plot_layout(heights = c(4,4,1.5)) +
    plot_annotation(theme = theme(plot.margin = unit(c(0.5, 3, 0, 0), "cm")))


ggsave("DNAkernels.png", plot = plt, path = "./figures/",
       width = 24, height = 20, units = "cm")



################################################################################

# Here using RNA to plot the NMDS plots for λ=0,1 for lienar, exponential and log
ncomps <- 4
pca.RNA <- pca(RNA.1 + 0.001, ncomp = ncomps, logratio = 'CLR')
RNAcomps <- pca.RNA$variates$X


# λ = 0
lambda <- 0
NoTime <- penaltyfunc(df = RNAcomps,
                      metric = CompBC,
                      reactor = reactors,
                      timepts = timepts,
                      kernelfunc = linearkernel,
                      lambda = lambda)
title <- expression(bold("No Reactor Time-Scale"))

MDSstructure <- mdsplotfunc(NoTime , conditions, timepts, title)
stress <- paste(100*round(MDSstructure$stress, 4), "%", sep = "")

MDS_NoTime  <- MDSstructure$mdsplot
MDS_NoTime  <- MDS_NoTime  + 
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    annotate("text", x = -Inf, y = -Inf, 
             label = stress, hjust = -0.15, vjust = -0.8,
             size = 5)


# Linear, λ = 1
lambda <- 1
Linear1 <- penaltyfunc(df = RNAcomps,
                       metric = CompBC,
                       reactor = reactors,
                       timepts = timepts,
                       kernelfunc = linearkernel,
                       lambda = lambda)
title <- "Linear Kernel"

MDSstructure <- mdsplotfunc(Linear1, conditions, timepts, title)
stress <- paste(100*round(MDSstructure$stress, 4), "%", sep = "")

MDS_Linear1 <- MDSstructure$mdsplot
MDS_Linear1 <- MDS_Linear1 + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    annotate("text", x = -Inf, y = -Inf, 
             label = stress, hjust = -0.15, vjust = -0.8,
             size = 5)



# Exponential, λ = 1
lambda <- 1
Exp1 <- penaltyfunc(df = RNAcomps,
                    metric = CompBC,
                    reactor = reactors,
                    timepts = timepts,
                    kernelfunc = expkernel,
                    lambda = lambda)
title <- "Exponential Kernel"

MDSstructure <- mdsplotfunc(Exp1, conditions, timepts, title)
stress <- paste(100*round(MDSstructure$stress, 4), "%", sep = "")

MDS_Exp1 <- MDSstructure$mdsplot
MDS_Exp1 <- MDS_Exp1 + 
    theme(legend.position = "none") +
    annotate("text", x = -Inf, y = -Inf, 
             label = stress, hjust = -0.15, vjust = -0.8,
             size = 5)


# Logarithmic, λ = 1
lambda <- 1
Log1 <- penaltyfunc(df = RNAcomps,
                    metric = CompBC,
                    reactor = reactors,
                    timepts = timepts,
                    kernelfunc = logkernel,
                    lambda = lambda)
title <- "Logarithmic Kernel"

MDSstructure <- mdsplotfunc(Log1, conditions, timepts, title)
stress <- paste(100*round(MDSstructure$stress, 4), "%", sep = "")

MDS_Log1 <- MDSstructure$mdsplot
MDS_Log1 <- MDS_Log1 + 
    theme(legend.position = c(1.15,1.1),
          axis.title.y = element_blank()) +
    annotate("text", x = -Inf, y = -Inf, 
             label = stress, hjust = -0.15, vjust = -0.8,
             size = 5)


plt <- (MDS_NoTime | MDS_Linear1) /
    (MDS_Exp1 | MDS_Log1) /
    wrapped_legend +
    plot_layout(heights = c(4,4, 1.5)) +
    plot_annotation(theme = theme(plot.margin = unit(c(0.5, 3, 0, 0), "cm")))


ggsave("RNAkernels.png", plot = plt, path = "./figures/",
       width = 24, height = 20, units = "cm")


################################################################################















