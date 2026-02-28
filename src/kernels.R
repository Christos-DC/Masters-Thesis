
# This is for the kernel functions implemented 

# Linear Kernel
linearkernel <- function(x, times){
    #' Inputs:
    #'      x: Non-negative value representing time difference
    #'      times: Vector of time points.
    #'      
    #' Returns: Linear Kernel value
    
    tstar <- max(times) - min(times)
    return(1 - x/tstar)
}

# Exponential Kernel
expkernel <- function(x, times) {
    #' Inputs:
    #'      x: Non-negative value representing time difference
    #'      times: Vector of time points.
    #'      
    #' Returns: Exponential value
    
    return(exp(-x))
    }


# Logarithmic Kernel
logkernel <- function(x, times){
    #' Inputs:
    #'      x: Non-negative value representing time difference
    #'      times: Vector of time points.
    #'      
    #' Returns: Logarithmic Kernel value
    
    tstar <- max(times) - min(times)
    return(1 - log(x+1)/log(tstar+1))
}



# Kernel matrix of the adjusted time points.
kernel_matrix <- function(kernelfunc, times){
    #' Inputs:
    #'      kernelfunc: The kernel function used for time adjustment
    #'      time: A vector of time points in sequential order
    #'      
    #' Returns: 
    #'      time_matrix: matrix of the kernel time lengths.
    
    # Setting up the time matrix
    n <- length(times)
    time_matrix <- matrix(0, nrow = n, ncol = n)
    diag(time_matrix) <- kernelfunc(0, times)
    
    # Gathering the kernel time point differences.
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            time_matrix[i,j] <- kernelfunc(abs(times[i] - times[j]), times)
            time_matrix[j,i] <- kernelfunc(abs(times[i] - times[j]), times)
        }
    }
    
    return(time_matrix)
}



# reactor: DNA3$reactor, timpts: DNA.3$week_from_salt1

# Adjusting the distance matrix based on the kernel function and lambda parameter.
penaltyfunc <- function(df, metric, reactor, timepts, kernelfunc, lambda){
    #' Inputs:
    #'      df: Raw abundance data frame
    #'      metric: distance metric to invoke on df.
    #'      reactor: Vector of reactors from df
    #'      timepts: Vector of time points from df
    #'      kernelfunc: The kernel function chosen for time adjustment.
    #'      lambda: Penalty value on the kernel function
    #'      
    #' Returns: Distance Matrix with adjusted kernel time-scaling.
    
    # Setting up parameters
    uniquereactors <- unique(reactor)
    n <- length(uniquereactors)
    
    # Gather the original distance metric
    Distdf <- sampledist(df, metric)
    Distdf <- data.frame(Distdf, reactor = reactor, time = timepts)
    
    index <- 0
    for (i in 1:n){
        # Filtering the distance matrix based on reactor
        react <- uniquereactors[i]
        subdf <- Distdf %>% filter(reactor == react)
        nsamp <- dim(subdf)[1]
        times <- subdf$time
        subdf <- subdf[colnames(subdf) %in% rownames(subdf), ]
        
        # Adjusting Distance values for reactor group.
        timemat <- kernel_matrix(kernelfunc, times)
        steps <- index + 1:nsamp
        subdist <- Distdf[steps, steps] * (1 - lambda * timemat)
        
        # Update the Distance Matrix
        Distdf[steps, steps] <- subdist
        index <- index + nsamp
    }
    
    return(Distdf)
}
