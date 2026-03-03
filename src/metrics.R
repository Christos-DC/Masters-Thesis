
#' This is all of the distance metrics that have been calculated.
#' 
#' Here the distance metrics are calculated where 0 and 1 means that the two
#' sample are completely the same and different respectively.
#' 
#' There is also a function outlining how I calculated the distance matrix.


# Calculates the Bray-Curtis value between two sample groups
Bray_Curtis <- function(sample1, sample2){
    #' Inputs:
    #'      sample1: A vector of integers representing the first sample group
    #'      sample2: A vector of integers representing the second sample group
    #' 
    #' Returns: 
    #'      Bray-Curtis value between two sample groups from 0 to 1 
    
    # Error: Plot lengths are different
    if (length(sample1) != length(sample2)){
        stop("Plots are not the same length.")
    }
    
    # Calculating the numerator and denominator terms
    numerator <- sum(abs(sample1 - sample2))
    denominator <- sum(sample1 + sample2)
    
    if (denominator == 0) {
        stop("Denominator is zero, cannot calculate Bray-Curtis dissimilarity.")
    }
    
    BC <- numerator / denominator
    return(BC)
}


# Calculates the PCA Bray-Curtis value between two sample groups.
PCA_Bray_Curtis <- function(sample1, sample2){
    #' Inputs:
    #'      sample1: A vector of integers representing the first sample group
    #'      sample2: A vector of integers representing the second sample group
    #' 
    #' Returns: 
    #'      PCA Bray-Curtis value between two sample groups from 0 to 1 
    
    # Error: Plot lengths are different
    if (length(sample1) != length(sample2)){
        stop("Plots are not the same length.")
    }
    
    # Calculating the numerator and denominator terms
    numerator <- sum(abs(sample1 - sample2))
    denominator <- sum(abs(sample1) + abs(sample2))
    
    if (denominator == 0) {
        stop("Denominator is zero, cannot calculate Bray-Curtis dissimilarity.")
    }
    
    BC <- numerator / denominator
    return(BC)
}


# Calculate the Jaccard value between two sample groups.
Jaccard <- function(sample1, sample2) {
    #' Inputs:
    #'      sample1: A vector of integers representing the first sample group
    #'      sample2: A vector of integers representing the second sample group
    #' 
    #' Returns: 
    #'      Jaccard value between two sample groups
    
    # Error: Plot lengths are different
    if (length(sample1) != length(sample2)){
        stop("Plots are not the same length.")
    }
    
    # Convert to binary and sum up the intersection and union values.
    intersection=sum(sample1>0 & sample2>0)
    union = sum(sample1>0) + sum(sample2>0) - intersection
    return (1-intersection/union)
}

# Calculate the Sorhensen-Dice value between two sample groups.
SDCoeff <- function(sample1, sample2){
    #' Inputs:
    #'      sample1: A vector of integers representing the first sample group
    #'      sample2: A vector of integers representing the second sample group
    #' 
    #' Returns: 
    #'      Sorhensen-Dice value between two sample groups from 0 to 1 
    
    # Error: Plot lengths are different
    if (length(sample1) != length(sample2)){
        stop("Plots are not the same length.")
    }
    
    # Convert to binary and sum up the intersection
    intersection=sum(sample1>0 & sample2>0)
    U <- sum(sample1 > 0);
    V <- sum(sample2 > 0);
    SD <- (2*intersection)/(U + V)
    
    return(1 - SD)
}


# Calculate the Normalised Expected Shared Species value between two sample groups.
NESS <- function(sample1, sample2, m = 5){
    #' Inputs:
    #'      sample1: A vector of integers representing the first sample group
    #'      sample2: A vector of integers representing the second sample group
    #'      m: positive integer number of random individuals drawn between the 
    #'         two groups.
    #' 
    #' Returns: 
    #'      Normalised Expected Shared Species value between two sample groups
    #'      from 0 to 1 
    
    # Error: m is not a positive integer
    if (m <= 0 || !is.integer(m)){
        stop("m is not a positive integer.")
    }
    
    # Error: Plot lengths are different
    if (length(sample1) != length(sample2)){
        stop("Plots are not the same length.")
    }
    
    # Calculating the Expected Shared Species
    ESS <- function(sample1, sample2, m){
        n1 <- sum(sample1)
        n2 <- sum(sample2)
        sp1 <- choose(n1 - sample1, m)
        sp2 <- choose(n2 - sample2, m)
        tot1 <- choose(n1,m)
        tot2 <- choose(n2,m)
        ESS_S1S2 <- sum((1- sp1/tot1) * (1 - sp2/tot2))
        return(ESS_S1S2)
    }
    
    # Calculate the Expected shared species for each group combination
    ESS_11 <- ESS(sample1, sample1, m)
    ESS_22 <- ESS(sample2, sample2, m)
    ESS_12 <- ESS(sample1, sample2, m)
    
    NESS = (2*ESS_12)/(ESS_11 + ESS_22)
    return(1-NESS)
}



# Hill Based Metric (Measures dissimilarity)
# Inspired by Oskar Modin's (2020) paper.
hill <- function(sample1, sample2, q=1){
    #' Inputs:
    #'      sample1: A vector of integers representing the first sample group
    #'      sample2: A vector of integers representing the second sample group
    #'      q: non-negative value that outlines the diversity order
    #' 
    #' Returns: 
    #'      Hill-based value between two sample groups from 0 to 1 
    
    # Error: q value chosen is less than zero.
    if (q < 0){
        stop("q cannot be less than zero")
    }
    
    # Error: Plot lengths are different.
    if (length(sample1) != length(sample2)){
        stop("Plots are not the same length.")
    }
    
    # scaling sample 1
    n1 <- sum(sample1)
    p1 <- sample1/n1
    
    # scaling sample 2
    n2 <- sum(sample2)
    p2 <- sample2/n2
    
    # gamma component for diversity calculation
    gam <- 0.5*p1 + 0.5*p2
    gam <- gam[gam > 0]
    p1 <- p1[p1 > 0]
    p2 <- p2[p2 > 0]
    
    if (q == 1){
        # Gamma Hill Number
        Hillgamma <- sum((gam) * log(gam))
        Hillgamma <- exp(-Hillgamma)
        
        # Alpha Hill Number
        Hillalpha <- 0.5*sum(p1 * log(p1)) + 0.5*sum(p2 * log(p2))
        Hillalpha <- exp(-Hillalpha)
        
        # Beta Hill number and metric
        Hillbeta <- Hillgamma / Hillalpha
        Hillmetric <- log(Hillbeta) / log(2)
    } 
    else {
        # Gamma Hill Number
        Hillgamma <- sum((gam)^q)
        Hillgamma <- Hillgamma^(1/(1-q))
        
        # Alpha Hill Number
        Hillalpha <- 0.5*sum(p1^q) + 0.5*sum(p2^q)
        Hillalpha <- Hillalpha^(1/(1-q))
        
        # Beta Hill number and metric
        Hillbeta <- Hillgamma / Hillalpha
        Hillmetric <- (Hillbeta^(1-q) - 1)/(2^(1-q) - 1)
    }
    
    return(Hillmetric)
}



# Calculating the distance between two samples depending on the metric function
sampledist <- function(X, dis.func, precision = 5, ...){
    #' Inputs:
    #'      X: is either a matrix or data frame that is the species abundance 
    #'         table; samples are in the rows and species are in the columns.
    #'      dis.func: function that calculates similarity.
    #'      precision: non-negative integer that outlines the decimal precision 
    #'                 of the final similarity values.
    #'      ... : other parameters that need to be passed for the dis.func. 
    #'            Example is the q parameter for the hill function.
    #' Returns:
    #'      Data frame of the Distance Matrix for all sample combinations.
    
    # Setting up the Distance Matrix
    n <- nrow(X)
    dis_table <- matrix(0, nrow = n, ncol = n)
    df_names <- rownames(X)
    
    # Calculate all possible dissimilarity combination values
    for (i in 1:(n-1)){ 
        for (j in (i+1):n){
            value <- dis.func(X[i,], X[j,], ...)
            value <- round(value, precision)
            dis_table[i, j] <- value
            dis_table[j, i] <- value
        }
    }
    
    # Renaming data frame
    dis_table <- as.data.frame(dis_table)
    rownames(dis_table) <- df_names
    colnames(dis_table) <- df_names
    
    return(dis_table)
}




#' # Testing a new metric
#' expMan <- function(sample1, sample2){
#'     #' Inputs:
#'     #'      sample1: A vector of integers representing the first sample group
#'     #'      sample2: A vector of integers representing the second sample group
#'     #'
#'     #' Returns:
#'     #'      Exponential Manhattan value between two sample groups
#' 
#'     # Error: Plot lengths are different
#'     if (length(sample1) != length(sample2)){
#'         stop("Plots are not the same length.")
#'     }
#' 
#'     Manhattan <- sum(abs(sample1 - sample2))
#'     SampleSum <- sum(abs(sample1) + abs(sample2))
#' 
#'     # Samples are completely distinct
#'     if (Manhattan == SampleSum) return(1)
#'     else return(1 - exp(-Manhattan))
#' }




# Second testing metric












