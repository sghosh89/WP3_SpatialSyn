rm(list = ls())

library(here)
library(tidyverse)

#===========================================================
# Bootstrap summaries: median, CI, p_boot
#===========================================================

boot_summarize_all <- function(fitboot, probs = c(0.025, 0.975)) {
  
  coefs <- stats::coef(fitboot)
  
  boot_mat <- fitboot$bootstrap
  
  if (is.null(boot_mat) || !is.matrix(boot_mat)) {
    stop("No bootstrap matrix found at fitboot$bootstrap.")
  }
  
  keep <- intersect(names(coefs), colnames(boot_mat))
  
  if (!length(keep)) {
    stop("No overlapping names between coef(fitboot) and bootstrap columns.")
  }
  
  boot_mat <- boot_mat[, keep, drop = FALSE]
  
  out <- lapply(keep, function(term) {
    
    draws <- boot_mat[, term]
    draws <- draws[is.finite(draws)]
    
    if (!length(draws)) {
      stop("All bootstrap draws are NA/Inf for term: ", term)
    }
    
    med <- stats::median(draws, na.rm = TRUE)
    ci  <- stats::quantile(draws, probs = probs, na.rm = TRUE)
    
    p_boot <- 2 * min(mean(draws >= 0), mean(draws <= 0))
    p_boot <- min(p_boot, 1)
    
    p_boot_se <- sqrt(p_boot * (1 - p_boot) / length(draws))
    
    data.frame(
      term      = term,
      slope_med = med,
      slope_lo  = unname(ci[1]),
      slope_hi  = unname(ci[2]),
      p_boot    = p_boot,
      p_boot_se = p_boot_se,
      estimate  = unname(coefs[term]),
      row.names = NULL
    )
  })
  
  do.call(rbind, out)
}

#===========================================================
# Add Wald SE, t-value, and p-value
#===========================================================

add_wald_pvals <- function(fitboot, boot_tab) {
  
  beta <- stats::coef(fitboot)
  V    <- if (!is.null(fitboot$vcov)) fitboot$vcov else stats::vcov(fitboot)
  se   <- sqrt(diag(V))
  
  n <- if (!is.null(fitboot$n)) fitboot$n else length(fitboot$y)
  d <- if (!is.null(fitboot$d)) fitboot$d else length(beta)
  
  df_res <- n - d
  
  terms <- intersect(boot_tab$term, names(beta))
  
  tval <- beta[terms] / se[terms]
  pval <- 2 * stats::pt(-abs(tval), df = df_res)
  
  wald_df <- data.frame(
    term      = terms,
    std_error = as.numeric(se[terms]),
    t_value   = as.numeric(tval),
    p_wald    = as.numeric(pval),
    row.names = NULL
  )
  
  out <- merge(
    boot_tab,
    wald_df,
    by = "term",
    all.x = TRUE,
    sort = FALSE
  )
  
  out <- out[match(boot_tab$term, out$term), , drop = FALSE]
  
  attr(out, "n")      <- n
  attr(out, "d")      <- d
  attr(out, "df_res") <- df_res
  
  out
}

#===========================================================
# Adaptive formatting: avoids 0.00 and -0.00
#===========================================================

pretty_num <- function(x,
                       digits = 3,
                       sci_digits = 2,
                       sci_threshold = 1e-3) {
  
  sapply(x, function(z) {
    
    if (is.na(z))
      return(NA_character_)
    
    if (z == 0)
      return("0")
    
    if (abs(z) < sci_threshold) {
      ## scientific notation
      return(formatC(z,
                     format = "e",
                     digits = sci_digits))
    } else {
      ## fixed decimal
      return(formatC(z,
                     format = "f",
                     digits = digits))
    }
    
  })
}

format_boot_table <- function(tab) {
  
  tab %>%
    mutate(
      across(
        c(
          slope_med,
          slope_lo,
          slope_hi,
          p_boot,
          p_boot_se,
          estimate,
          std_error,
          t_value,
          p_wald
        ),
        pretty_num
      )
    )
}

#===========================================================
# Helper function: run one model and save summary
#===========================================================

run_boot_summary <- function(resloc) {
  
  fitboot <- readRDS(
    paste0(resloc, "/model_est_phylolm_boot.RDS")
  )
  
  boot_tab  <- boot_summarize_all(fitboot)
  boot_tab2 <- add_wald_pvals(fitboot, boot_tab)
  boot_tab3 <- format_boot_table(boot_tab2)
  
  sink(
    here(paste(resloc, "boot_summary.txt", sep = "/")),
    append = FALSE,
    split = TRUE
  )
  
  print(boot_tab3)
  
  cat("\n")
  cat("n      =", attr(boot_tab2, "n"), "\n")
  cat("d      =", attr(boot_tab2, "d"), "\n")
  cat("df_res =", attr(boot_tab2, "df_res"), "\n")
  cat("\n")
  
  sink()
  
  invisible(boot_tab3)
}

#===========================================================
# Model folders
#===========================================================

model_dirs <- c(
  
  # 0-250 km, 32 years
  "RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5_abs_td",
  "RESULTS/model_phylolm_sig95_0-250km_min32yr/model_tas5_abs_td",
  "RESULTS/model_phylolm_sig95_0-250km_min32yr/model_tas5",
  "RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5",
  
  # 0-250 km, 36 years
  "RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5_abs_td",
  "RESULTS/model_phylolm_sig95_0-250km_min36yr/model_tas5_abs_td",
  "RESULTS/model_phylolm_sig95_0-250km_min36yr/model_tas5",
  "RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5",
  
  # 0-250 km, 40 years
  "RESULTS/model_phylolm_sig95_0-250km/model_pr5_abs_td",
  "RESULTS/model_phylolm_sig95_0-250km/model_tas5_abs_td",
  "RESULTS/model_phylolm_sig95_0-250km/model_tas5",
  "RESULTS/model_phylolm_sig95_0-250km/model_pr5",
  
  # 0-100 km, 40 years
  "RESULTS/model_phylolm_sig95_0-100km/model_pr5_abs_td",
  "RESULTS/model_phylolm_sig95_0-100km/model_tas5_abs_td",
  
  # 100-250 km, 40 years
  "RESULTS/model_phylolm_sig95_100-250km/model_pr5_abs_td",
  "RESULTS/model_phylolm_sig95_100-250km/model_tas5_abs_td",
  "RESULTS/model_phylolm_sig95_100-250km/model_tas5",
  "RESULTS/model_phylolm_sig95_100-250km/model_pr5"
)

#===========================================================
# Run all summaries
#===========================================================

all_boot_summaries <- lapply(model_dirs, function(x) {
  run_boot_summary(here(x))
})

names(all_boot_summaries) <- model_dirs

#================= 
source(here("R/phylolm_model_with_sampling_effort.R"))
#===================

#===========================================================
# Sampling-effort sensitivity models
#===========================================================

sampling_model_dirs <- c(
  
  # 32 years
  "RESULTS/model_phylolm_sig95_0-250km_min32yr",
  
  # 36 years
  "RESULTS/model_phylolm_sig95_0-250km_min36yr",
  
  # 40 years
  "RESULTS/model_phylolm_sig95_0-250km"
  
)

run_boot_sampling_summary <- function(resloc, climate="pr"){
  
  fitbootP <- readRDS(paste0(resloc, "/model_est_phylolm_boot_with_samplingeffort_P.RDS"))
  fitbootT <- readRDS(paste0(resloc, "/model_est_phylolm_boot_with_samplingeffort_T.RDS"))
  
  boot_tabP  <- boot_summarize_all(fitbootP)
  boot_tabT  <- boot_summarize_all(fitbootT)
  
  sink(
    here(paste(resloc, "boot_summary_samplingeffort.txt", sep = "/")),
    append = FALSE,
    split = TRUE
  )
  
  print(boot_tabP)
  
  cat("\n")
  cat("\n")
  cat("\n")
  
  print(boot_tabT)
  
  sink()
}


sampling_boot_summaries <- lapply(sampling_model_dirs, function(x) {
  run_boot_sampling_summary(here(x))
})

names(sampling_boot_summaries) <- sampling_model_dirs


