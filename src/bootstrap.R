
# This is for the Bootstrapping components in the Thesis

euclidean <- function(x,y){
    #' Inputs:
    #'      x, y: n-length vector of real numbers.
    #'      
    #'  Returns: Euclidean distance between x and y.
    
    if (length(x) != length(y)){
        stop("Vectors are of different length!")
    }
    return(sqrt(sum((x-y)^2)))
}


BlockBootstrap <- function(df, b){
    #' Inputs:
    #'      df: A data frame with ordered structure along the rows. 
    #'      b: Positive integer that indicates the block size.
    #'      
    #' Returns: A block bootstrapped data frame from df.
    
    df <- as.data.frame(df)
    
    # Define the parameters
    n <- dim(df)[1]
    space <- floor(b/2)
    low <- space + 1
    high <- n - space 
    
    # Invalid choice of b
    if (low >= high){
        stop("b value is not suitable!")
    }
    
    # Randomise the ranges for a Bootstrapped dataset 
    indexrange <- low:high
    indices <- sample(indexrange, ceiling(n/b), replace = TRUE)
    newdf <- NULL
    
    # Odd case
    if (b %% 2 == 1){
        for (index in indices){
            adjustdf <- as.data.frame(df[(index - space):(index + space), ])
            colnames(adjustdf) <- colnames(df)
            newdf <- rbind(newdf, adjustdf)
        }
    }
    else{ # Even case
        for (index in indices){
            bernoulli <- rbinom(1,1,1/2)
            
            # Extra block to the right
            if (bernoulli == 1){
                adjustdf <- as.data.frame(df[(index - space + 1):(index + space), ])
                colnames(adjustdf) <- colnames(df)
                newdf <- rbind(newdf, adjustdf)
            } # Extra block to the left
            else{
                adjustdf <- as.data.frame(df[(index - space):(index + space - 1), ])
                colnames(adjustdf) <- colnames(df)
                newdf <- rbind(newdf, adjustdf)
            }
        }
    }
    
    # Truncate the final values to keep the data size the same.
    newdf <- newdf[1:n, ] 
    
    return(newdf)
}

MDSnumeric <- function(df, func, reactor, condition, timepts, distdf = FALSE,  
                       boot = TRUE, ...){
    #' Inputs:
    #'      df: A data frame for the original abundance data or distance matrix
    #'      func: A ecological metric of choice.
    #'      reactor: vector of the reactors each sample groups are in.
    #'      condition: vector of the condition each sample groups are in.
    #'      timepts: vector of the time points each sample groups belong too.
    #'      distdf: FALSE if df is raw abundance data. TRUE if df is distance matrix
    #'      boot: Boolean value for performing blocked bootstrapping.
    #'      ... : Additional parameters passed for the metric chosen.
    #' 
    #' Returns: A List containing the following
    #'      meanvals: Vector of the mean values for each condition
    #'      sdvals (boot = TRUE): Vector of the bootstrapped standard deviations for each condition.
    #'      medianvals: Vector of the median values for each condition
    #'      confint (boot = TRUE): Data frame of the 95% bootstrapped confidence intervals.
    #'      Nconditions: Number of data points for each condition.
    #'      time: The time it took to run the function.
    
    # Start the time 
    start <- Sys.time()
    
    # Get the mds data frame
    if (distdf){
        MDSdf <- mds(df, type = "ordinal", ndim = 2)
    }
    else {
        distmat <- sampledist(df, func, ...)
        MDSdf <- mds(distmat, type = "ordinal", ndim = 2)
    }
    dfReduction <- MDSdf$conf
    dfReduction <- data.frame(dfReduction, "Reactor" = reactor, 
                              "Condition" = condition, "Time" = timepts)
    
    calcEuclidean <- function(df, n){
        #' Inputs:
        #'      df: Data frame of the coordinates from mds result.
        #'      n: A positive integer of the number of points in df
        #' 
        #' Returns: A vector of weighted average distance values.
        
        Distvals <- c()
        
        # Calculating weighted time average
        for (i in 2:(n-1)){
            # Distance between current point and its previous and future point.
            distvaldown <- euclidean(df[i, c(1,2)], df[i-1, c(1,2)])
            distvalup <- euclidean(df[i+1, c(1,2)], df[i, c(1,2)])
            
            # Gathering time weights
            timetotal <- abs(timepts[i+1] - timepts[i-1])
            timedown <- abs(timepts[i] - timepts[i-1])
            timeup <- abs(timepts[i+1] - timepts[i])
            
            # Calculating weighted average.
            distval <- timedown * distvaldown + timeup * distvalup
            Distvals <- c(Distvals, distval)
        }
        
        return(Distvals)
    }
    
    
    # Store the Statistics
    UniqueConditions <- sort(unique(condition))
    Nconditions <- length(UniqueConditions)
    NdistValues <- means <- sds <- medians <- rep(0, Nconditions) 
    ConfIntervals <- matrix(0, nrow = Nconditions, ncol = 2)
    
    # How many bootstrap samples we want to include.
    N <- 1000 
    for (k in 1:Nconditions){
        # Filter data based on the condition
        condition_df <- dfReduction %>% filter(Condition == UniqueConditions[k])
        UniqueReactors <- unique(condition_df$Reactor)
        Nreactors <- length(UniqueReactors)
        
        Distvals <- NULL
        simDistvalsCond <- NULL
        
        for (j in 1:Nreactors){
            # Further filtering for each reactor in one condition group.
            reactor_df <- condition_df %>% filter(Reactor == UniqueReactors[j])
            reactor_df <- reactor_df %>% arrange(Time)
            n <- dim(reactor_df)[1]
            NdistValues[k] <- NdistValues[k] + n - 2
            
            # Calculating the distance metric
            reactDistvals <- calcEuclidean(reactor_df[,c(1,2)], n)
            Distvals <- c(Distvals, reactDistvals)
            
            # If bootstrapping is set to TRUE calculate bootstrap samples
            if (boot){
                simDistvalsReact <- NULL 
                for (i in 1:N){
                    # Choosing block size of 2 
                    simVals<- BlockBootstrap(reactDistvals, 2) 
                    simVals <- t(as.data.frame(simVals))
                    simDistvalsReact <- rbind(simDistvalsReact, simVals)
                }
                simDistvalsCond <- cbind(simDistvalsCond, simDistvalsReact)
            }
        }
        
        # Collect statistics
        means[k] <- mean(Distvals)
        medians[k] <- median(Distvals)
        
        # Gather bootstrapped statistics
        if (boot){
            avgDistvals <- rowMeans(simDistvalsCond)
            sds[k] <- sd(avgDistvals)
            ConfIntervals[k, ] <- quantile(avgDistvals, c(0.025, 0.975))
        }
    }
    
    
    # Naming data structures
    names(NdistValues) <- UniqueConditions
    names(means) <- UniqueConditions
    names(sds) <- UniqueConditions
    names(medians) <- UniqueConditions
    colnames(ConfIntervals) <- c("2.5%", "97.5%")
    rownames(ConfIntervals) <- UniqueConditions
    
    end <- Sys.time() - start
    # Return the appropriate list based on if bootstrapped is performed.
    if (boot){
        return(list(meanvals = means, sdvals = sds, medianvals = medians,
                    confint = ConfIntervals, Nconditions = NdistValues,
                    time = end))
    }
    else{
        return(list(meanvals = means, medianvals = medians, 
                    Nconditions = NdistValues, time = end))
    }
}









