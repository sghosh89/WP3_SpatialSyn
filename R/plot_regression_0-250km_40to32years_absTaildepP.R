#==================== plot across species ==============
rm(list=ls())

library(here)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(readxl)

nbin <- 4

#==================== read data ====================

dfsig95.40 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km/df_sub.RDS"))
dfsig95.32 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min32yr/df_sub.RDS"))
dfsig95.36 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min36yr/df_sub.RDS"))

dfsig95.40.abs <- dfsig95.40 %>%
  dplyr::select(abs.avg.td.pr5.sig, fpr5.sig,
                abs.avg.td.tas5.sig, ftas5.sig,
                abs.avg.td.ab.sig, fab.sig, HWI) %>%
  mutate(year = "40 years")

dfsig95.36.abs <- dfsig95.36 %>%
  dplyr::select(abs.avg.td.pr5.sig, fpr5.sig,
                abs.avg.td.tas5.sig, ftas5.sig,
                abs.avg.td.ab.sig, fab.sig, HWI) %>%
  mutate(year = "36 years")

dfsig95.32.abs <- dfsig95.32 %>%
  dplyr::select(abs.avg.td.pr5.sig, fpr5.sig,
                abs.avg.td.tas5.sig, ftas5.sig,
                abs.avg.td.ab.sig, fab.sig, HWI) %>%
  mutate(year = "32 years")

dfsig.abs <- rbind(dfsig95.40.abs, dfsig95.36.abs, dfsig95.32.abs)

levels_year <- c("40 years", "36 years", "32 years")

dfsig.abs$year <- factor(dfsig.abs$year, levels = levels_year)

pal <- c(
  "40 years" = "#000000",
  "36 years" = "#cc79a7",
  "32 years" = "#A6761D"
)

#==================== read bootstrap draws: precipitation model ====================

boot32 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5_abs_td/model_est_phylolm_boot.RDS"))
boot36 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5_abs_td/model_est_phylolm_boot.RDS"))
boot40 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km/model_pr5_abs_td/model_est_phylolm_boot.RDS"))

x32 <- as.data.frame(boot32$bootstrap) %>%
  mutate(year = "32 years")

x36 <- as.data.frame(boot36$bootstrap) %>%
  mutate(year = "36 years")

x40 <- as.data.frame(boot40$bootstrap) %>%
  mutate(year = "40 years")

xall <- rbind(x40, x36, x32)

xall$year <- factor(xall$year, levels = levels_year)

source(here("R/helper_plotter.R"))
#==================== precipitation panel ====================

g1 <- ggplot(
  data = dfsig.abs,
  aes(x = abs.avg.td.pr5.sig,
      y = abs.avg.td.ab.sig,
      col = year)
) +
  geom_point(pch = 19, alpha = 0.5) +
  xlab("Average (absolute) tail-dependent spatial synchrony,\n Precipitation") +
  ylab("Average (absolute) tail-dependent spatial synchrony,\n Abundance") +
  theme_bw(base_size = 20) +
  facet_wrap(~year, scale = "free") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank()
  )

# Adjusted prediction:
# precipitation varies, HWI fixed at year-specific mean
pred_df_P <- make_adjusted_pred(
  data = dfsig.abs,
  boot_draws = xall,
  xvar = "abs.avg.td.pr5.sig",
  covar = "HWI"
)

g1P <- g1 +
  geom_line(
    data = pred_df_P,
    aes(x = x, y = y50, color = year),
    size = 1,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df_P,
    aes(x = x, y = ylo, color = year),
    linetype = "dashed",
    size = 1,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df_P,
    aes(x = x, y = yhi, color = year),
    linetype = "dashed",
    size = 1,
    inherit.aes = FALSE
  )

label_P <- make_label_df(
  file40 = "RESULTS/model_phylolm_sig95_0-250km/model_pr5_abs_td/boot_summary.txt",
  file36 = "RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5_abs_td/boot_summary.txt",
  file32 = "RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5_abs_td/boot_summary.txt",
  focal_term = "abs.avg.td.pr5.sig"
)

g1P <- g1P +
  geom_label(
    data = label_P,
    aes(
      x = -Inf,
      y = Inf,
      label = label,
      color = year
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 5,
    fill = NA,
    label.size = 0.25,
    show.legend = FALSE
  )

#==================== HWI panel ====================

g1.hwi <- ggplot(
  data = dfsig.abs,
  aes(x = HWI,
      y = abs.avg.td.ab.sig,
      col = year)
) +
  geom_point(pch = 19, alpha = 0.5) +
  xlab("Hand-Wing Index, HWI") +
  ylab("Average (absolute) tail-dependent spatial synchrony,\n Abundance") +
  theme_bw(base_size = 20) +
  facet_wrap(~year, scale = "free") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank()
  )

# Adjusted prediction:
# HWI varies, precipitation fixed at year-specific mean
pred_df_HWI <- make_adjusted_pred(
  data = dfsig.abs,
  boot_draws = xall,
  xvar = "HWI",
  covar = "abs.avg.td.pr5.sig"
)

g1P.hwi <- g1.hwi +
  geom_line(
    data = pred_df_HWI,
    aes(x = x, y = y50, color = year),
    size = 1,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df_HWI,
    aes(x = x, y = ylo, color = year),
    linetype = "dashed",
    size = 1,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df_HWI,
    aes(x = x, y = yhi, color = year),
    linetype = "dashed",
    size = 1,
    inherit.aes = FALSE
  )

label_HWI <- make_label_df(
  file40 = "RESULTS/model_phylolm_sig95_0-250km/model_pr5_abs_td/boot_summary.txt",
  file36 = "RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5_abs_td/boot_summary.txt",
  file32 = "RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5_abs_td/boot_summary.txt",
  focal_term = "HWI"
)

g1P.hwi <- g1P.hwi +
  geom_label(data = label_HWI,
             aes(x = -Inf, y = Inf, label = label, color = year),
             inherit.aes = FALSE,
             hjust = -0.05,
             vjust = 1.1,
             size = 5,
             fill = NA,
             label.size = 0.25,
             show.legend = FALSE)

#==================== save figure ====================

pdf(here("Results/plot_regression_0-250km_40to32years_absTaildep_P.pdf"), width = 18, height = 12)
grid.arrange(g1P, g1P.hwi)
dev.off()

# caption: 
#   Figure X. Relationships between average absolute tail-dependent spatial synchrony in abundance and 
# average absolute tail-dependent spatial synchrony in precipitation (top row) and Hand–Wing Index (HWI; bottom row) 
# across species under three minimum time-series length thresholds (40, 36, and 32 years). 
# Points represent species-level observations. Solid lines show the median adjusted fitted relationship estimated 
# from phylogenetic linear models (phylolm, λ model) using 1000 parametric bootstrap replicates. 
# For each panel, predictions were generated by varying the focal predictor while holding the non-focal covariate 
# at its year-specific sample mean (i.e., precipitation effects were evaluated at mean HWI and HWI effects at mean precipitation synchrony). 
# Dashed lines indicate the 95% bootstrap interval of the fitted relationship, 
# calculated as the 2.5th and 97.5th quantiles of predicted mean responses across bootstrap coefficient draws. 
# Insets report the bootstrap median slope estimate, bootstrap 95% confidence interval of the focal regression coefficient, 
# and bootstrap p-value from the multivariable phylogenetic regression model:  Abundance synchrony ~ precipitation synchrony + HWI.






