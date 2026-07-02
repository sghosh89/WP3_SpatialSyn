# before proceeding with phylolm, or pgls, try a linear fit model just to see if there is any outlier
# this could arise because of species which were rarely sampled
# if you screened such species, maybe you could put a certain threshold that because of precision, (lowering variance/ or uncertainty),
# we will be confident about out tail-dep. estimate if species are sampled, say, at least 10 sites (i.e., nint = 45 unique site pairs)

# initial diagonistics

# Are species-level tail estimates (abs.avg.td.ab.sig, abs.avg.td.pr5.sig, abs.avg.td.tas5.sig) dependent on nint?
library(here)
library(ggplot2)
library(patchwork)
library(dplyr)

plot_nint_sensitivity <- function(yr_threshold,
                                    site_thr = 10,
                                    sitepair_thr = 45){
  
  ##------------------------------------------------------------
  ## Read data
  ##------------------------------------------------------------
  outputresloc <- if (yr_threshold == 40){
    here("RESULTS/model_phylolm_sig95_0-250km")
  }else{
    here(paste0("RESULTS/model_phylolm_sig95_0-250km_min",
                yr_threshold,"yr"))
  }
  
  df <- readRDS(file.path(outputresloc,"df.RDS"))
  
  tmp <- df %>%
    filter(!is.na(abs.avg.td.ab.sig),
           !is.na(nint),
           nint > 0) %>%
    mutate(
      abs_dev = abs(
        abs.avg.td.ab.sig -
          mean(abs.avg.td.ab.sig, na.rm = TRUE)
      )
    )
  
  ##------------------------------------------------------------
  ## Axis breaks
  ##------------------------------------------------------------
  
  ## site numbers to display on the top axis
  site_breaks <- c(2,3,4,5,6,8,10,15,20,30,50,100)
  
  ## convert sites -> unique site pairs
  pair_breaks <- site_breaks*(site_breaks-1)/2
  
  ##------------------------------------------------------------
  ## Common theme
  ##------------------------------------------------------------
  
  common_theme <-
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title.x.top = element_text(margin = margin(b = 8))
    )
  
  ##------------------------------------------------------------
  ## Tail estimate
  ##------------------------------------------------------------
  
  p1 <-
    ggplot(tmp,
           aes(x = nint,
               y = abs.avg.td.ab.sig)) +
    
    geom_point(size = 2, alpha = 0.8) +
    
    geom_smooth(method = "lm",
                se = TRUE,
                colour = "blue") +
    
    geom_vline(xintercept = sitepair_thr,
               linetype = "dashed") +
    
    annotate("text",
             x = sitepair_thr,
             y = Inf,
             label = "10 sites",
             angle = 90,
             vjust = -0.4,
             hjust = 1.1,
             fontface = "bold",
             size = 3.4) +
    
    scale_x_log10(
      breaks = pair_breaks,
      labels = pair_breaks,
      sec.axis = dup_axis(
        breaks = pair_breaks,
        labels = site_breaks,
        name = "Sampled sites"
      )
    ) +
    
    common_theme +
    
    labs(
      title = paste0("sampled for a minimum of ", yr_threshold, " years"),
      x = "Unique site pairs (log scale)",
      y = expression(A[TD]^abundance)
    )
  
  ##------------------------------------------------------------
  ## Absolute deviation
  ##------------------------------------------------------------
  
  p2 <-
    ggplot(tmp,
           aes(x = nint,
               y = abs_dev)) +
    
    geom_point(size = 2, alpha = 0.8) +
    
    geom_smooth(method = "loess",
                se = TRUE,
                colour = "red") +
    
    geom_vline(xintercept = sitepair_thr,
               linetype = "dashed") +
    
    annotate("text",
             x = sitepair_thr,
             y = Inf,
             label = "10 sites",
             angle = 90,
             vjust = -0.4,
             hjust = 1.1,
             fontface = "bold",
             size = 3.4) +
    
    scale_x_log10(
      breaks = pair_breaks,
      labels = pair_breaks,
      sec.axis = dup_axis(
        breaks = pair_breaks,
        labels = site_breaks,
        name = "Sampled sites"
      )
    ) +
    
    common_theme +
    
    labs(
      title = "Absolute deviation from mean",
      x = "Unique site pairs (log scale)",
      y = "Absolute deviation"
    )
  
  final_plot <- p1 / p2
  
  ggsave(
    file.path(outputresloc,
              "nint_sensitivity_check.pdf"),
    final_plot,
    width = 6.5,
    height = 7
  )
  
  final_plot
}

## Run
p40 <- plot_nint_sensitivity(40)
p36 <- plot_nint_sensitivity(36)
p32 <- plot_nint_sensitivity(32)


# =========== for 40 years: plot interpretation ===================
#1. Mean dependence (top row)
# 
# Goal:
# 
# Does the average tail estimate change with nint?
# 
# Abundance
# 
# There is a strong negative trend, but visually it is almost entirely driven by one extreme point at very low nint (your species with only 6 pairs).
# 
# Without that point, the cloud looks much flatter.
# 
# Precipitation
# 
# Looks approximately flat.
# I do not see strong evidence that mean tail dependence changes with nint.
# 
# Temperature
# 
# Maybe a weak positive trend, but uncertainty is large.
# 
# So from the top row alone:
# 
# No obvious systematic mean dependence after scaling
# Except possible influence of the low-nint abundance outlier
# 2. Variance dependence (bottom row)
# 
# Goal:
# 
# Are estimates less precise when nint is small?
# 
# This is the more important result.
# 
# Abundance
# 
# Very strong pattern:
# 
# enormous spread at low nint
# rapidly stabilizes
# 
# This is exactly what we expected:
# species with few site pairs produce noisy estimates.
# 
# Precipitation
# 
# Moderate heteroscedasticity.
# Variance appears somewhat larger at low nint.
# 
# Temperature
# 
# Similar story; uncertainty larger at low nint.
# 
# So:
# 
# Var(
# L
# ^
# )↓asnint↑
# 
#================= for 36 years : plot interpretation ==================
# 1. Mean dependence (top row)
# 
# Question:
#   
#   Does average tail dependence still depend on nint after dividing by nint?
#   
#   Abundance
# 
# There is still a negative trend visually, but it looks driven by a few low-nint observations. Most species are concentrated and nearly flat at larger nint.
# 
# Precipitation
# 
# Looks essentially flat.
# 
# Temperature
# 
# Very weak positive tendency, but small.
# 
# Overall:
#   
#   I do not see strong evidence that mean tail estimates systematically increase with nint.
# 
# That supports your normalization.
# 
# 2. Variance dependence (bottom row)
# 
# Question:
#   
#   Are estimates less precise at low nint?
#   
#   Abundance
# 
# Very strong.
# 
# Variance collapses rapidly as nint increases.
# 
# This is exactly what you would expect if estimates based on few site pairs are noisy.
# 
# Precipitation
# 
# Moderate reduction in spread.
# 
# Temperature
# 
# Similar but weaker.
# 
# So the message becomes:
#   
#   Tail estimates appear approximately unbiased with respect to nint, but uncertainty decreases substantially as the number of site pairs increases.
# 
#================= comparing with 32 years =========
#Comparing thresholds
#Threshold	Pattern
#40 yr	One low-nint species dominates
#36 yr	Much cleaner
#32 yr	Cleanest, most stable

#That does not necessarily mean 32 years is better biologically, but statistically it suggests:

#more retained species
#more site pairs
#reduced influence of sparse estimates

# To evaluate whether species-level tail-dependence estimates were sensitive to differences in spatial sampling, 
# we examined relationships between normalized tail estimates and the number of unique site pairs (nint). 
# Mean tail estimates showed little dependence on nint, indicating that normalization by the number of site pairs 
# reduced sampling-related bias. However, variability in estimated tail dependence declined with increasing nint, 
# suggesting lower precision for species represented by few site pairs.




