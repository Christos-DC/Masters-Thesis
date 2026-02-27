
# This is for the Bootstrapping components in the Thesis

euclidean <- function(x,y){
    if (length(x) != length(y)){
        stop("Vectors are of different length!")
    }
    return(sqrt(sum((x-y)^2)))
}


BlockBootstrap <- function(df, b){
    # df should be a data frame of temporal data points via rows. So t times r data frame.
    # b is the block size of the bootstrapping.
    df <- as.data.frame(df)
    
    n <- dim(df)[1]
    space <- floor(b/2)
    low <- space + 1
    high <- n - space 
    
    if (low >= high){
        stop("chosen b is not suitable!")
    }
    
    indexrange <- low:high
    indices <- sample(indexrange, ceiling(n/b), replace = TRUE)
    newdf <- NULL
    if (b %% 2 == 1){
        for (index in indices){
            adjustdf <- as.data.frame(df[(index - space):(index + space), ])
            colnames(adjustdf) <- colnames(df)
            newdf <- rbind(newdf, adjustdf)
        }
    }
    else{ # Splitting the even cases depending on whether accepting blocked values going up or down
        for (index in indices){
            bernoulli <- rbinom(1,1,1/2)
            if (bernoulli == 1){
                adjustdf <- as.data.frame(df[(index - space + 1):(index + space), ])
                colnames(adjustdf) <- colnames(df)
                newdf <- rbind(newdf, adjustdf)
            }
            else{
                adjustdf <- as.data.frame(df[(index - space):(index + space - 1), ])
                colnames(adjustdf) <- colnames(df)
                newdf <- rbind(newdf, adjustdf)
            }
        }
    }
    newdf <- newdf[1:n, ] # truncate the final values to keep the data size the same (not perfect but is okay for now).
    
    return(newdf)
}

MDSnumeric <- function(df, func, reactor, condition, timepts, distdf = FALSE, timeweight = FALSE, 
                       boot = TRUE,...){
    # df is the original data frame of interest.
    # func is the distance metric used to calculate the dissimilarity.
    # reactor is about what reactor the sample groups are in.
    # condition is about the type of group each reactor belongs too.
    # timepts are the time points for the experiment.
    # distdf is used if the df is an actual distance matrix or the raw abundance data.
    # timeweight is a true/false variable that outlines whether the results account for the time distance between points or not.
    # ... for any other parameters needed for the sampledist function.
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
        Distvals <- c()
        # Weighted time average
        if (timeweight){
            for (i in 2:(n-1)){
                distvaldown <- euclidean(df[i, c(1,2)], df[i-1, c(1,2)])
                distvalup <- euclidean(df[i+1, c(1,2)], df[i, c(1,2)])
                
                timetotal <- abs(timepts[i+1] - timepts[i-1])
                timedown <- abs(timepts[i] - timepts[i-1])
                timeup <- abs(timepts[i+1] - timepts[i])
                
                distval <- timedown * distvaldown + timeup * distvalup
                Distvals <- c(Distvals, distval)
            }
        } # Standard time average. (not sure if I need to keep this since I want to reflect the actual time differences. Probably can ask Saritha about this since I think it's redundant, but am keeping this for now to see if this properly works.)
        else{
            for (i in 2:(n-1)){
                distvaldown <- euclidean(df[i, c(1,2)], df[i-1, c(1,2)])
                distvalup <- euclidean(df[i+1, c(1,2)], df[i, c(1,2)])
                distval <- (distvaldown + distvalup)/2
                Distvals <- c(Distvals, distval)
            }
        }
        return(Distvals)
    }
    
    
    # Store the Statistics
    UniqueConditions <- sort(unique(condition))
    Nconditions <- length(UniqueConditions)
    NdistValues <- means <- sds <- medians <- rep(0, Nconditions) 
    ConfIntervals <- matrix(0, nrow = Nconditions, ncol = 2)
    
    
    N <- 1000 # How many bootstrap samples we want to include.
    for (k in 1:Nconditions){
        condition_df <- dfReduction %>% filter(Condition == UniqueConditions[k])
        UniqueReactors <- unique(condition_df$Reactor)
        Nreactors <- length(UniqueReactors)
        
        Distvals <- NULL
        simDistvalsCond <- NULL
        # Further filtering for each reactor in one condition group.
        for (j in 1:Nreactors){
            reactor_df <- condition_df %>% filter(Reactor == UniqueReactors[j])
            reactor_df <- reactor_df %>% arrange(Time)
            n <- dim(reactor_df)[1]
            NdistValues[k] <- NdistValues[k] + n - 2
            
            # Calculating the distance metric
            reactDistvals <- calcEuclidean(reactor_df[,c(1,2)], n)
            Distvals <- c(Distvals, reactDistvals)
            
            if (boot){
                simDistvalsReact <- NULL 
                for (i in 1:N){
                    simVals<- BlockBootstrap(reactDistvals, 2) #choosing 2 because of the way the distance is calculated.
                    simVals <- t(as.data.frame(simVals))
                    simDistvalsReact <- rbind(simDistvalsReact, simVals)
                }
                simDistvalsCond <- cbind(simDistvalsCond, simDistvalsReact)
            }
        }
        
        means[k] <- mean(Distvals)
        medians[k] <- median(Distvals)
        
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