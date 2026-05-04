
# Here are the functions that are used to plot the MDS functions

# Function that blends the colours over time.
blend_with_white <- function(colours, timepts){
    #' Inputs:
    #'      colours: vector that has colour inputs.
    #'      timepts: vector that consists of time points (not unique timepts)
    #'      
    #' Returns:
    #'      The blended rgb colours.
    
    n <- length(timepts)
    
    # Preparing the rgb values and the white matrix
    rgb_vals <- matrix(0, nrow = n, ncol = 3)
    for (i in 1:n){
        rgb_vals[i,] <- t(col2rgb(colours[i]) / 255)
    }
    white <- matrix(1, nrow = n, ncol = 3)
    
    # If negative time points, shift the scale down by min(temp) - 1
    if (min(timepts) <= 0){
        timepts <- timepts - min(timepts) +1
    }
    
    # Fading the colours based on the time points
    max_timept <- max(timepts)
    ratio <- timepts / max_timept
    blended <- ratio * rgb_vals + (1-ratio) * white
    
    # Convert to HTML colour code
    rgb_cols <- rep(NA, n)
    for (i in 1:n){
        rgb_cols[i] <- rgb(blended[i,1], blended[i,2], blended[i,3])
    }
    
    return(rgb_cols)
}

# The base colours for the three groups
base_colours <- c("control" = "red", "inhib_1" = "blue", "inhib_2" = "green")

# MDS plot function for the Zymo dataset
mdsplotfunc <- function(df, condition, timepts, title){
    #' Inputs:
    #'      df: Data Frame of the distance matrix.
    #'      condition: vector outlines the conditions the data points are in.
    #'      timepts: The time points of when the abundance data was collected.
    #'      title: String that outlines the title of the plot.
    #'      
    #' Return: A list containing
    #'      mdsplot: The MDS plot of the distance matrix.
    #'      stress: A positive float value outlining the stress value.
    
    # Perform MDS on the distance matrix
    mdsresult <- mds(df, type = "ordinal")
    mdsresult_conf <- as.data.frame(cbind(mdsresult$conf, "time" = timepts, "condition" = condition))
    mdsresult_conf[,1:3] <- sapply(mdsresult_conf[,1:3], as.numeric)
    
    # Gather the colours for the different groups
    hex_colors <- blend_with_white(base_colours[mdsresult_conf$condition], mdsresult_conf$time)
    
    # Plotting the results
    mdsplot <- mdsresult_conf %>% 
        mutate(color = hex_colors) %>% 
        ggplot(aes(x=D1, y=D2, shape = condition)) +
        geom_point(aes(color = color), size = 3) +
        labs(title = title, shape = "Condition") +
        scale_color_identity() +
        theme(legend.title = element_text("Condition")) +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA),
              panel.border = element_rect(fill = NA, color = "black"),
              plot.title = element_text(size = 12, hjust = 0.5))
    
    # Getting the stress of the mds plot
    mdsStress <- mdsresult$stress
    
    return(list(mdsplot=mdsplot, stress=mdsStress))
}








