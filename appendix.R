
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
    df <- penaltyfunc(comps, PCA_Bray_Curtis, reactors, timepts, linearkernel, lambda)
    DNAnumerics <- MDSnumeric(df = df[, -c(99,100)], 
                              func = PCA_Bray_Curtis, 
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
    df <- penaltyfunc(comps, PCA_Bray_Curtis, reactors, timepts, linearkernel, lambda)
    DNAnumerics <- MDSnumeric(df = df[, -c(99,100)], 
                              func = PCA_Bray_Curtis, 
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





































