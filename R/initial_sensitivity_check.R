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

plot_nint_sensitivity <- function(yr_threshold, sitepair_thr=45){
  
  if(yr_threshold==40){
    outputresloc<-here("RESULTS/model_phylolm_sig95_0-250km")
  }else{
    outputresloc<-here(paste0("RESULTS/model_phylolm_sig95_0-250km_min",yr_threshold,"yr"))
  }
  df<-readRDS(here(paste(outputresloc,"/df.RDS",sep="")))
  
  vars <- c(
    "abs.avg.td.ab.sig",
    "abs.avg.td.pr5.sig",
    "abs.avg.td.tas5.sig"
  )
  
  titles <- c(
    "Abundance",
    "Precipitation",
    "Temperature"
  )
  
  p_mean <- list()
  p_var <- list()
  
  for(i in seq_along(vars)){
    
    v <- vars[i]
    
    tmp <- df %>%
      filter(
        !is.na(.data[[v]]),
        !is.na(nint),
        nint > 0
      )
    
    # ------------------------
    # Row 1 → Mean check
    # ------------------------
    
    p_mean[[i]] <-
      ggplot(
        tmp,
        aes(
          x = log10(nint),
          y = .data[[v]]
        )
      ) +
      geom_point(size=2) +
      geom_smooth(
        method="lm",
        se=TRUE,
        color="blue"
      ) +
      geom_vline(
        xintercept = log10(sitepair_thr),
        linetype = "dashed"
      )+
      #annotate("text", x = log10(sitepair_thr),
       # y = Inf, label = paste0("n=", sitepair_thr),
       # vjust = 1.5, angle = 90)+
      theme_bw() +
      labs(
        title=paste("Mean:", titles[i]),
        x="log10(unique site pairs)",
        y="Tail estimate"
      )
    
    
    # ------------------------
    # Row 2 → Variance check
    # ------------------------
    
    tmp$dev <-
      abs(
        tmp[[v]] -
          mean(
            tmp[[v]],
            na.rm=TRUE
          )
      )
    
    p_var[[i]] <-
      ggplot(
        tmp,
        aes(
          x=log10(nint),
          y=dev
        )
      ) +
      geom_point(size=2) +
      geom_smooth(
        method="loess",
        se=TRUE,
        color="red"
      ) +
      geom_vline(
        xintercept = log10(sitepair_thr),
        linetype = "dashed"
      )+
      theme_bw() +
      labs(
        title=paste("Variance:", titles[i]),
        x="log10(unique site pairs)",
        y="Absolute deviation"
      )
    
  }
  
  final_plot <-
    (p_mean[[1]] | p_mean[[2]] | p_mean[[3]]) /
    (p_var[[1]]  | p_var[[2]]  | p_var[[3]])
  
  ggsave(here(paste(outputresloc,"/nint_sensitivity_check.pdf",sep="")),
         final_plot,
         width=12,
         height=7
  )
  
  return(final_plot)
  
}

yr_threshold<-40
p40 <- plot_nint_sensitivity(yr_threshold=yr_threshold)


yr_threshold<-36
p36 <- plot_nint_sensitivity(yr_threshold=yr_threshold)

yr_threshold<-32
p32 <- plot_nint_sensitivity(yr_threshold=yr_threshold)


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

#==================
library(car)
detect_outliers_cook <- function(yr_threshold, formula, species_col = "Species"){
  
  if(yr_threshold==40){
    outputresloc<-here("RESULTS/model_phylolm_sig95_0-250km")
  }else{
    outputresloc<-here(paste0("RESULTS/model_phylolm_sig95_0-250km_min",yr_threshold,"yr"))
  }
  df<-readRDS(here(paste(outputresloc,"/df.RDS",sep="")))
  
  fit <- lm(formula, data = df)
  
  cook <- cooks.distance(fit)
  
  thr <- 4 / nrow(df)
  
  res <- data.frame(Species = df[[species_col]],
    cooks = cook,
    influential = cook > thr
  )
  
  # plot
  plot(cook, pch = 19,
       ylab = "Cook's distance", 
       xlab = "Observation",
       main = deparse(formula))
  
  abline(h = thr, col = "red", lty = 2, lwd = 2)
  
  text(which(cook > thr),
    cook[cook > thr],
    labels = df[[species_col]][cook > thr],
    pos = 3,
    cex = 0.7
  )
  
  cat("\nCook threshold =", round(thr,3), "\n\n")
  
  print(
    res[order(-res$cooks),]
  )
  
  return(res)
  
}

yr_threshold<-40
res_ab40 <- detect_outliers_cook(yr_threshold=yr_threshold, formula= abs.avg.td.ab.sig ~ log10(nint), species_col = "Species")

yr_threshold<-36
res_ab36 <- detect_outliers_cook(yr_threshold=yr_threshold, formula= abs.avg.td.ab.sig ~ log10(nint), species_col = "Species")

yr_threshold<-32
res_ab32 <- detect_outliers_cook(yr_threshold=yr_threshold, formula= abs.avg.td.ab.sig ~ log10(nint), species_col = "Species")

#====================================================
compare_fit_without_outlier <- function(yr_threshold, formula){
  
  if(yr_threshold==40){
    outputresloc<-here("RESULTS/model_phylolm_sig95_0-250km")
  }else{
    outputresloc<-here(paste0("RESULTS/model_phylolm_sig95_0-250km_min",yr_threshold,"yr"))
  }
  df<-readRDS(here(paste(outputresloc,"/df.RDS",sep="")))
  df<-df%>%filter(nint>=45)
  
  # Fit full model
  fit <- lm(formula, data=df)
  
  print(summary(fit))
}

yr_threshold<-40
res40 <- compare_fit_without_outlier(yr_threshold=yr_threshold,  abs.avg.td.ab.sig ~ abs.avg.td.pr5.sig)
resf40 <- compare_fit_without_outlier(yr_threshold=yr_threshold,  fab.sig ~ fpr5.sig)

