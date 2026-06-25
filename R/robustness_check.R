#===========================================================
# Leave-one-out diagnostics for phylolm regression models
# Models:
#   Abundance synchrony ~ Precipitation synchrony + HWI
#   Abundance synchrony ~ Temperature synchrony + HWI
# Across:
#   40 years, 36 years, 32 years
#===========================================================

rm(list = ls())

library(tidyverse)
library(phylolm)
library(ape)
library(here)
library(gridExtra)

pal_pred <- c(
  "Precipitation" = "#0072B2",   # blue
  "Temperature"   = "#A62800"    # red
)

#===========================================================
# General leave-one-out function for phylolm
#===========================================================

run_loo_phylolm <- function(data, phy, formula, focal_term,
                            model = "lambda",
                            year_label,
                            predictor_label) {
  
  phy$node.label <- NULL
  
  if (!"Species" %in% names(data)) {
    stop("The data must contain a column named 'Species'.")
  }
  
  rownames(data) <- data$Species
  
  spp_keep <- intersect(rownames(data), phy$tip.label)
  data <- data[spp_keep, , drop = FALSE]
  phy <- drop.tip(phy, setdiff(phy$tip.label, spp_keep))
  data <- data[phy$tip.label, , drop = FALSE]
  
  n <- nrow(data)
  
  fit_full <- phylolm(
    formula = formula,
    data = data,
    phy = phy,
    model = model
  )
  
  full_sum <- summary(fit_full)
  full_slope <- coef(fit_full)[focal_term]
  full_p <- full_sum$coefficients[focal_term, "p.value"]
  
  loo_res <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    removed_species <- rownames(data)[i]
    data_i <- data[-i, , drop = FALSE]
    phy_i <- drop.tip(phy, removed_species)
    
    fit_i <- tryCatch(
      phylolm(
        formula = formula,
        data = data_i,
        phy = phy_i,
        model = model
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_i)) {
      
      loo_res[[i]] <- data.frame(
        year = year_label,
        predictor = predictor_label,
        removed = removed_species,
        slope = NA_real_,
        p_value = NA_real_,
        lambda = NA_real_
      )
      
    } else {
      
      fit_sum <- summary(fit_i)
      
      loo_res[[i]] <- data.frame(
        year = year_label,
        predictor = predictor_label,
        removed = removed_species,
        slope = coef(fit_i)[focal_term],
        p_value = fit_sum$coefficients[focal_term, "p.value"],
        lambda = fit_i$optpar
      )
    }
  }
  
  loo_res <- bind_rows(loo_res)
  
  loo_summary <- loo_res %>%
    summarise(
      year = year_label,
      predictor = predictor_label,
      full_slope = full_slope,
      full_p = full_p,
      loo_min_slope = min(slope, na.rm = TRUE),
      loo_median_slope = median(slope, na.rm = TRUE),
      loo_max_slope = max(slope, na.rm = TRUE),
      prop_same_sign = mean(sign(slope) == sign(full_slope), na.rm = TRUE),
      prop_p_lt_0.05 = mean(p_value < 0.05, na.rm = TRUE),
      prop_p_lt_0.10 = mean(p_value < 0.10, na.rm = TRUE),
      most_influential_removed = removed[which.max(abs(slope - full_slope))]
    )
  
  return(
    list(
      full_model = fit_full,
      loo_results = loo_res,
      loo_summary = loo_summary
    )
  )
}

#===========================================================
# Load data and trees
#===========================================================

df40 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km/df_sub.RDS"))
df36 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min36yr/df_sub.RDS"))
df32 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min32yr/df_sub.RDS"))

ct3_40 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km/consensus_tree_with_edgelength.RDS"))
ct3_36 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min36yr/consensus_tree_with_edgelength.RDS"))
ct3_32 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min32yr/consensus_tree_with_edgelength.RDS"))

#===========================================================
# Run LOO models
#===========================================================

loo_T_40 <- run_loo_phylolm(
  data = df40,
  phy = ct3_40,
  formula = abs.avg.td.ab.sig ~ abs.avg.td.tas5.sig + HWI,
  focal_term = "abs.avg.td.tas5.sig",
  year_label = "40 years",
  predictor_label = "Temperature"
)

loo_T_36 <- run_loo_phylolm(
  data = df36,
  phy = ct3_36,
  formula = abs.avg.td.ab.sig ~ abs.avg.td.tas5.sig + HWI,
  focal_term = "abs.avg.td.tas5.sig",
  year_label = "36 years",
  predictor_label = "Temperature"
)

loo_T_32 <- run_loo_phylolm(
  data = df32,
  phy = ct3_32,
  formula = abs.avg.td.ab.sig ~ abs.avg.td.tas5.sig + HWI,
  focal_term = "abs.avg.td.tas5.sig",
  year_label = "32 years",
  predictor_label = "Temperature"
)

loo_P_40 <- run_loo_phylolm(
  data = df40,
  phy = ct3_40,
  formula = abs.avg.td.ab.sig ~ abs.avg.td.pr5.sig + HWI,
  focal_term = "abs.avg.td.pr5.sig",
  year_label = "40 years",
  predictor_label = "Precipitation"
)

loo_P_36 <- run_loo_phylolm(
  data = df36,
  phy = ct3_36,
  formula = abs.avg.td.ab.sig ~ abs.avg.td.pr5.sig + HWI,
  focal_term = "abs.avg.td.pr5.sig",
  year_label = "36 years",
  predictor_label = "Precipitation"
)

loo_P_32 <- run_loo_phylolm(
  data = df32,
  phy = ct3_32,
  formula = abs.avg.td.ab.sig ~ abs.avg.td.pr5.sig + HWI,
  focal_term = "abs.avg.td.pr5.sig",
  year_label = "32 years",
  predictor_label = "Precipitation"
)

#===========================================================
# Combine results
#===========================================================

loo_all_results <- bind_rows(
  loo_T_40$loo_results,
  loo_T_36$loo_results,
  loo_T_32$loo_results,
  loo_P_40$loo_results,
  loo_P_36$loo_results,
  loo_P_32$loo_results
) %>%
  mutate(
    year = factor(year, levels = c("40 years", "36 years", "32 years")),
    predictor = factor(predictor, levels = c("Precipitation", "Temperature")),
    significant = p_value < 0.05
  )

loo_all_summary <- bind_rows(
  loo_T_40$loo_summary,
  loo_T_36$loo_summary,
  loo_T_32$loo_summary,
  loo_P_40$loo_summary,
  loo_P_36$loo_summary,
  loo_P_32$loo_summary
) %>%
  mutate(
    year = factor(year, levels = c("40 years", "36 years", "32 years")),
    predictor = factor(predictor, levels = c("Precipitation", "Temperature"))
  )

#===========================================================
# Save tables
#===========================================================

write.csv(
  loo_all_results,
  here("RESULTS/LOO_phylolm_all_results.csv"),
  row.names = FALSE
)

write.csv(
  loo_all_summary,
  here("RESULTS/LOO_phylolm_summary.csv"),
  row.names = FALSE
)

#===========================================================
# Summary table for p-value robustness
#===========================================================

loo_p_summary <- loo_all_results %>%
  group_by(predictor, year) %>%
  summarise(
    n_loo = sum(!is.na(p_value)),
    median_p = median(p_value, na.rm = TRUE),
    min_p = min(p_value, na.rm = TRUE),
    max_p = max(p_value, na.rm = TRUE),
    prop_p_lt_0.05 = mean(p_value < 0.05, na.rm = TRUE),
    prop_p_lt_0.10 = mean(p_value < 0.10, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  loo_p_summary,
  here("RESULTS/LOO_phylolm_pvalue_summary.csv"),
  row.names = FALSE
)

print(loo_all_summary)
print(loo_p_summary) # see the prop_p_lt_0.05 column; show this in suppmat table

#===========================================================
# Diagnostic plot 1: LOO slope stability
#===========================================================

p_slope <- ggplot(
  loo_all_results,
  aes(x = year, y = slope, color = predictor)
) +
  geom_jitter(width = 0.15, alpha = 0.65, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ predictor, scales = "fixed") +
  scale_color_manual(values = pal_pred) +
  theme_bw(base_size = 16) +
  ylab("Leave-one-out slope estimate") +
  xlab("Minimum time-series length") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

#===========================================================
# Diagnostic plot 2: LOO p-value stability
#===========================================================

p_pvalue <- ggplot(
  loo_all_results,
  aes(x = year, y = p_value, color = predictor)
) +
  geom_jitter(
    aes(shape = significant),
    width = 0.15,
    alpha = 0.75,
    size = 2
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ predictor, scales = "fixed") +
  scale_color_manual(values = pal_pred) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "p < 0.05", `FALSE` = "p > 0.05")
  ) +
  theme_bw(base_size = 16) +
  ylab("Leave-one-out p-value") +
  xlab("Minimum time-series length") +
  labs(shape = "LOO significance") +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

#===========================================================
# Save diagnostic plots
#===========================================================

pdf(
  here("RESULTS/LOO_phylolm_diagnostic_slope_pvalue.pdf"),
  width = 16,
  height = 10
)

grid.arrange(
  p_slope,
  p_pvalue,
  ncol = 1
)

dev.off()

#===========================================================
# Show plots in RStudio
#===========================================================

print(p_slope)
print(p_pvalue)

#Leave-one-out sensitivity analysis showed that significance of precipitation effects was preserved across nearly all refits (>95%), 
#whereas temperature effects became robustly detectable only under the 32-year threshold.

#=========== methods ======
# Sensitivity analysis using leave-one-out phylogenetic regression
# 
# To evaluate whether regression estimates were disproportionately influenced by individual species, we performed a 
# leave-one-out (LOO) sensitivity analysis for each phylogenetic linear model (phylolm, λ model). For each fitted model, 
# one species was removed at a time, the corresponding phylogenetic tree was pruned, and the model was re-estimated using the remaining species. 
# This procedure was repeated for all species within each dataset and separately for each minimum time-series length threshold (40, 36, and 32 years). 
# For every refit, we recorded the estimated regression coefficient of the focal predictor and its associated p-value.
# 
# We summarized LOO results using diagnostic distributions of (i) leave-one-out slope estimates and (ii) leave-one-out p-values. 
# Slope diagnostics were used to evaluate directional stability of the estimated relationship, 
# whereas p-value diagnostics quantified robustness of statistical significance to species removal. 
# We additionally computed the proportion of LOO refits remaining significant (p < 0.05) as a summary measure of inference stability.

#============== results ==========

# Leave-one-out analyses showed contrasting robustness patterns between precipitation and temperature effects. 
# For precipitation synchrony, regression slopes remained consistently positive across all temporal thresholds, 
# and statistical significance was preserved in nearly all leave-one-out refits. This indicates that the observed 
# precipitation–abundance synchrony relationship was not driven by a small number of influential species.
# 
# For temperature synchrony, leave-one-out slopes remained directionally consistent across thresholds; however, 
# statistical significance emerged only under the 32-year minimum time-series criterion. Importantly, 
# the significance observed for the 32-year subset persisted across nearly all leave-one-out refits, 
# indicating that the detected temperature effect was not attributable to one or a few influential observations. 
# In contrast, the 40- and 36-year subsets remained predominantly non-significant under species 
# deletion despite maintaining similar slope directions, suggesting insufficient effect 
# strength or precision rather than outlier dependence.
# 
# Overall, these results indicate that precipitation effects are robust across temporal thresholds, 
# whereas temperature effects appear threshold-dependent and become reliably detectable only under the 32-year criterion.


