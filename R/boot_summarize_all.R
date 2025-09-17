rm(list=ls())
library(here)
library(tidyverse)

#' Bootstrap summaries (median/CI) + two-sided p_boot for all coefficients
#' Works for phylolm objects fitted with boot > 0

# Bootstrap summaries + two-sided p_boot for all regression coefficients
boot_summarize_all <- function(fitboot, probs = c(0.025, 0.975)) {
  # 1) get coefficient point estimates
  coefs <- stats::coef(fitboot)                  # named numeric vector
  
  # 2) get bootstrap matrix (your object has it as a numeric matrix)
  boot_mat <- fitboot$bootstrap
  if (is.null(boot_mat) || !is.matrix(boot_mat)) {
    stop("No bootstrap matrix found at fitboot$bootstrap.")
  }
  
  # 3) keep only columns that correspond to coefficients (drop sigma2, lambda)
  keep <- intersect(names(coefs), colnames(boot_mat))
  if (!length(keep)) stop("No overlapping names between coef(fitboot) and bootstrap columns.")
  boot_mat <- boot_mat[, keep, drop = FALSE]
  
  # 4) compute stats per coefficient
  out <- lapply(keep, function(term) {
    draws <- boot_mat[, term]
    draws <- draws[is.finite(draws)]
    if (!length(draws)) stop("All bootstrap draws are NA/Inf for term: ", term)
    
    med <- stats::median(draws, na.rm = TRUE)
    ci  <- stats::quantile(draws, probs = probs, na.rm = TRUE)
    
    # two-sided bootstrap p-value: mass on the opposite side of 0 (×2)
    p_boot <- 2 * min(mean(draws >= 0), mean(draws <= 0))
    
    # Monte Carlo SE for p_boot (useful when B is modest)
    mc_se <- sqrt(p_boot * (1 - p_boot) / length(draws))
    
    data.frame(
      term       = term,
      #draws      = length(draws),
      slope_med  = med,
      slope_lo  = unname(ci[1]),
      slope_hi   = unname(ci[2]),
      p_boot     = p_boot,
      p_boot_se  = mc_se,
      estimate   = unname(coefs[term]),
      row.names  = NULL
    )
  })
  
  do.call(rbind, out)
}

# Add Wald p-values from summary(fitboot) to an existing bootstrap table
# Add Wald p-values computed from coef() and vcov() to an existing bootstrap table
# Add Wald SE, t-value, df, and t-based p-value (matches summary(phylolm))
# Add Wald SE, t, and t-based p-value; also expose n, d, and df_res = n - d
add_wald_pvals <- function(fitboot, boot_tab) {
  beta <- stats::coef(fitboot)
  V    <- if (!is.null(fitboot$vcov)) fitboot$vcov else stats::vcov(fitboot)
  se   <- sqrt(diag(V))
  
  # compute n, d, and residual df the same way phylolm does
  n  <- if (!is.null(fitboot$n)) fitboot$n else length(fitboot$y)
  d  <- if (!is.null(fitboot$d)) fitboot$d else length(beta)
  df_res <- n - d
  
  # align with rows in boot_tab
  terms <- intersect(boot_tab$term, names(beta))
  tval  <- beta[terms] / se[terms]
  pval  <- 2 * stats::pt(-abs(tval), df = df_res)
  
  wald_df <- data.frame(
    term      = terms,
    std_error = as.numeric(se[terms]),
    t_value   = as.numeric(tval),
    #df_res    = df_res,
    p_wald    = as.numeric(pval),
    row.names = NULL
  )
  
  out <- merge(boot_tab, wald_df, by = "term", all.x = TRUE, sort = FALSE)
  out <- out[match(boot_tab$term, out$term), , drop = FALSE]
  
  # attach handy attributes so you can fetch them later
  attr(out, "n")      <- n
  attr(out, "d")      <- d
  attr(out, "df_res") <- df_res
  out
}


#==================================
# now calling the function:
# 0-250km distance; 32 year

resloc<-here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)


sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_tas5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)


sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_tas5")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()


#-----------------------------
# now calling the function:
# 0-250km distance; 36 year

resloc<-here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_tas5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_tas5")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

#-----------------------------
# now calling the function:
# 0-250km distance; 40 year

resloc<-here("RESULTS/model_phylolm_sig95_0-250km/model_pr5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km/model_tas5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km/model_tas5")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)


sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-250km/model_pr5")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

#-----------------------------
# now calling the function:
# 0-100km distance; 40 year

resloc<-here("RESULTS/model_phylolm_sig95_0-100km/model_pr5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_0-100km/model_tas5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

# now calling the function:
# 100-250km distance; 40 year

resloc<-here("RESULTS/model_phylolm_sig95_100-250km/model_pr5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_100-250km/model_tas5_abs_td")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_100-250km/model_tas5")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

resloc<-here("RESULTS/model_phylolm_sig95_100-250km/model_pr5")
fitboot32 <- readRDS(paste0(resloc,"/model_est_phylolm_boot.RDS"))
boot_tab <- boot_summarize_all(fitboot32)
boot_tab2 <- add_wald_pvals(fitboot32, boot_tab)

sink(here(paste(resloc,"boot_summary.txt",sep="/")),append=TRUE, split=TRUE)
print(boot_tab2)
attr(boot_tab2, "n"); attr(boot_tab2, "d"); attr(boot_tab2, "df_res")
sink()

