
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
base_colors <- c("control" = "red", "inhib_1" = "blue", "inhib_2" = "green")

# MDS plot function for the Zymo dataset
mdsplotfunc <- function(df, title){
    #' Inputs:
    #'      df: Data Frame of the distance matrix.
    #'      title: String that outlines the title of the plot.
    #'      
    #' Return: A list containing
    #'      mdsplot: The MDS plot of the distance matrix.
    #'      stress: A positive float value outlining the stress value.
    
    # Perform MDS on the distance matrix
    mdsresult <- mds(df, type = "ordinal")
    mdsresult_conf <- as.data.frame(cbind(mdsresult$conf, "time" = DNA.3$week_from_salt1, "condition" = DNA.3$condition))
    mdsresult_conf[,1:3] <- sapply(mdsresult_conf[,1:3], as.numeric)
    
    # Gather the colours for the different groups
    hex_colors <- blend_with_white(base_colors[mdsresult_conf$condition], mdsresult_conf$time)
    
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
              plot.title = element_text(size = 10))
    
    # Getting the stress of the mds plot
    mdsStress <- mds$stress
    
    return(list(mdsplot=mdsplot, stress=mdsStress))
}



# Plotting the legend plot for Control, Inhibition 1 and 2.
mintime <- min(DNA.3$week_from_salt1)
maxtime <- max(DNA.3$week_from_salt1)
n <- 100
weeks <- round(seq(mintime, maxtime, length.out = n), 4)

# Gather the blended colours over time.
legend_df <- expand.grid(sort(names(base_colors)), weeks)
colnames(legend_df) <- c("condition", "week")
legend_df$color <- blend_with_white(base_colors[legend_df$condition], as.numeric(legend_df$week)) 

# Legend plot (compact, horizontal)
legend_plot <- ggplot(legend_df, aes(x = week, y = condition, fill = color)) +
    geom_tile(height = 0.6) +
    scale_fill_identity() +
    theme_minimal(base_size = 7) +
    labs(x = "Week", y = NULL) +
    scale_x_continuous(breaks = seq(-3,11, by=1)) +
    theme(
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", color = NA)
    )

# Control the dimension of the plot
wrapped_legend <- wrap_elements(legend_plot) +
    plot_layout(widths = unit(0.5, "npc"))





















