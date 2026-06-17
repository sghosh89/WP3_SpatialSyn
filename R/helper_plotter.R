#==================== helper function ====================
# This creates adjusted bootstrap predictions from:
# A = Intercept + beta_x * x + beta_covar * mean(covar)

#===========================================================
# Function: make_adjusted_pred()
#
# Purpose:
#   Create adjusted fitted lines from a phylogenetic multiple regression:
#
#       A = Intercept + beta_x * X + beta_z * Z
#
#   where:
#     X = focal predictor shown on the x-axis
#     Z = non-focal covariate held fixed at its year-specific mean
#
#   Example 1:
#     For precipitation panel:
#       X = precipitation synchrony
#       Z = HWI
#
#   Example 2:
#     For HWI panel:
#       X = HWI
#       Z = precipitation synchrony
#
#   The function uses all bootstrap coefficient draws to calculate:
#     y50 = median fitted line
#     ylo = 2.5% bootstrap line
#     yhi = 97.5% bootstrap line

make_adjusted_pred <- function(data, boot_draws, xvar, covar) {
  
  #---------------------------------------------------------
  # 1. Get the observed range of the focal x-axis variable
  #    separately for each year threshold.
  #    This ensures the fitted line is drawn only across
  #    the observed data range in each facet.
  #---------------------------------------------------------
  
  x_ranges <- data %>%
    group_by(year) %>%
    summarise(
      x_min = min(.data[[xvar]], na.rm = TRUE),
      x_max = max(.data[[xvar]], na.rm = TRUE),
      .groups = "drop"
    )
  
  #---------------------------------------------------------
  # 2. Compute the year-specific mean of the non-focal
  #    covariate.
  #
  #    This is the key multiple-regression correction.
  #    Instead of ignoring the other predictor, we hold it
  #    fixed at its typical value within each year subset.
  #---------------------------------------------------------
  
  covar_means <- data %>%
    group_by(year) %>%
    summarise(
      covar_mean = mean(.data[[covar]], na.rm = TRUE),
      .groups = "drop"
    )
  
  #---------------------------------------------------------
  # 3. Use each bootstrap coefficient draw to generate
  #    predicted abundance synchrony across a grid of x values.
  #
  #    For each bootstrap draw b:
  #
  #      y_b = Intercept_b + beta_x,b * x + beta_z,b * mean(Z)
  #
  #    where:
  #      beta_x,b = bootstrap slope for focal predictor
  #      beta_z,b = bootstrap slope for non-focal covariate
  #---------------------------------------------------------
  
  pred_df <- boot_draws %>%
    
    # Keep only year, intercept, focal coefficient,
    # and non-focal coefficient from bootstrap output
    dplyr::select(year, `(Intercept)`, all_of(xvar), all_of(covar)) %>%
    
    # Create a bootstrap draw ID within each year
    group_by(year) %>%
    mutate(draw = row_number()) %>%
    ungroup() %>%
    
    # Attach observed x-axis range for that year
    inner_join(x_ranges, by = "year") %>%
    
    # Attach year-specific mean of non-focal covariate
    inner_join(covar_means, by = "year") %>%
    
    # For each year, create 200 evenly spaced x values
    # between observed min and max
    group_by(year) %>%
    group_modify(~ {
      crossing(
        .x,
        x = seq(.x$x_min[1], .x$x_max[1], length.out = 200)
      )
    }) %>%
    ungroup() %>%
    
    # Multiple-regression-adjusted prediction:
    #
    # predicted abundance synchrony =
    #   intercept
    # + focal slope * focal x value
    # + non-focal slope * mean(non-focal covariate)
    mutate(
      y = `(Intercept)` +
        .data[[xvar]] * x +
        .data[[covar]] * covar_mean
    ) %>%
    
    # At each year and x value, summarize the 1000 bootstrap
    # fitted values into median and 95% bootstrap interval
    group_by(year, x) %>%
    summarise(
      y50 = median(y, na.rm = TRUE),
      ylo = quantile(y, 0.025, na.rm = TRUE),
      yhi = quantile(y, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(pred_df)
}

#==================== helper function for labels ====================

make_label_df <- function(file40, file36, file32, focal_term) {
  
  tab40 <- read.table(here(file40), header = TRUE, nrows = 3)
  tab36 <- read.table(here(file36), header = TRUE, nrows = 3)
  tab32 <- read.table(here(file32), header = TRUE, nrows = 3)
  
  s40 <- tab40[tab40$term == focal_term, ]
  s36 <- tab36[tab36$term == focal_term, ]
  s32 <- tab32[tab32$term == focal_term, ]
  
  tab_all <- rbind(s40, s36, s32)
  tab_all <- as.data.frame(tab_all)
  tab_all$year <- factor(levels_year, levels = levels_year)
  
  tab_all <- tab_all %>%
    mutate(
      label = sprintf(
        "Slope Median = %.2f\nSlope 95%% CI: [%.2f, %.2f]\np_boot = %.3f",
        slope_med, slope_lo, slope_hi, p_boot
      )
    )
  
  return(tab_all)
}
