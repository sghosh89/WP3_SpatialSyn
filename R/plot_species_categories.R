library(ggplot2)

plot_species_categories <- function(
    yr_threshold = 40,
    title = "Species Counts Across Categories A-E",
    palette = "Set2"
) {
  
  # Select counts based on threshold
  counts <- switch(
    as.character(yr_threshold),
    "40" = c(173, 95, 78, 35, 34),
    "36" = c(276, 106, 170, 93, 86),# note: 85 sp used in phylogenetic analysis as from BirdTree name matching
    "32" = c(328, 88, 240, 136, 124),# note: 122 sp used in phylogenetic analysis as from BirdTree name matching
    stop("yr_threshold must be 32, 36, or 40")
  )
  
  # Category labels with descriptions
  labels <- c(
    "A\nSampled ≥2 sites",
    "B\nNo spatial\nsynchrony",
    "C\nSpatial\nsynchrony",
    "D\nSpatial synchrony\nwith significant\ntail-dependence",
    "E\nCategory D\nsampled at\n≥10 sites"
  )
  
  # Create dataframe
  species_data <- data.frame(
    Category = factor(c("A", "B", "C", "D", "E"),
                      levels = c("A", "B", "C", "D", "E")),
    Number_of_species = counts
  )
  
  # Plot
  p <- ggplot(
    species_data,
    aes(
      x = Category,
      y = Number_of_species,
      fill = Category
    )
  ) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = Number_of_species),
      vjust = -0.5,
      size = 5
    ) +
    scale_x_discrete(labels = labels) +
    scale_fill_brewer(palette = palette) +
    labs(
      x = NULL,
      y = "Number of species",
      title = paste0(title, " (sites were sampled within 250 km for a minimum of ", yr_threshold, " years)")
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(
        size = 15,
        color = "black",
        lineheight = 0.9
      ),
      axis.text.y = element_text(
        size = 15,
        color = "black",
        lineheight = 0.9
      ),
      plot.title = element_text(
        hjust = 0.5,
        color = "black",
        face = "bold"
      )
    )
  
  return(p)
}

# Example
#plot_species_categories(yr_threshold = 40)
#plot_species_categories(yr_threshold = 36)
#plot_species_categories(yr_threshold = 32)