# --- Packages ---
rm(list=ls())
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

#--------------------

source(here("R/plot_species_categories.R"))


visualize_spat_syn_Nyrs_threshold <- function(yr_threshold = 40, target_dist_cat = c(0, 250), 
                                              nbin = 4, siglevel = "95%CI", label_size = 6){
  
  if (yr_threshold == 40) {
    df1 <- readRDS(here(paste(
      "RESULTS/abundance_spatsyn_nbin_", nbin,
      "_corlmcoru_sigres_summary_", target_dist_cat[1], "-",
      target_dist_cat[2], "Km.RDS", sep = ""
    )))
  } else {
    df1 <- readRDS(here(paste(
      "RESULTS/abundance_spatsyn_nbin_", nbin,
      "_corlmcoru_sigres_summary_", target_dist_cat[1], "-",
      target_dist_cat[2], "Km_min", yr_threshold, "yr.RDS", sep = ""
    )))
  }
  
  if (siglevel == "95%CI") {
    df1 <- df1 %>%
      filter(Lsig95ab != 0 | Usig95ab != 0) %>%
      dplyr::select(
        AOU,
        nint,
        Lsigab = Lsig95ab,
        Usigab = Usig95ab
      )
  }
  
  ## Compute percentage contributions
  df1 <- df1 %>%
    filter(nint >= 45) %>%
    mutate(
      lower_strength = Lsigab / nint,
      upper_strength = abs(Usigab) / nint,
      total_strength = lower_strength + upper_strength,
      fLval = 100 * lower_strength / total_strength,
      fUval = 100 * upper_strength / total_strength
    ) %>%
    arrange(desc(fLval))          # Red → Blue
  
  #nrow(df1)# 34 for 40 years, 86 for 36 years, 124 for 32 years
  #sum(df1$fLval>50.0)# 15 out of 34 for 40 yrs, 29 out of 86 for 36 years, 40 out of 124 for 32 years
  #sum(df1$fLval<50.0)# 19 out of 34 for 40 yrs, 57 out of 86 for 36 years, 84 out of 124 for 32 years
  
  #binom.test(15, 34, p = 0.5)
  #binom.test(29, 86, p = 0.5, alternative = "less")# p<0.001 means reject H0, accept HA
  #binom.test(40, 124, p = 0.5)
  
  ## Order species
  ord_vec <- df1$AOU
  
  ## Convert to long format
  df <- df1 %>%
    select(individual = AOU, fLval, fUval) %>%
    pivot_longer(
      cols = c(fLval, fUval),
      names_to = "observation",
      values_to = "value"
    ) %>%
    mutate(
      individual = factor(individual,
                          levels = rev(ord_vec)),   # red at top after coord_flip()
      observation = factor(
        observation,
        levels = c("fLval", "fUval"),
        labels = c("Lower-tail dependence",
                   "Upper-tail dependence")
      )
    )
  
  ## Plot
  p <- ggplot(df,
              aes(x = individual,
                  y = value,
                  fill = observation)) +
    geom_col(position = "fill",
             width = 0.8,
             alpha = 0.95) +
    coord_flip() +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1),
      labels = c(0, 25, 50, 75, 100),
      expand = c(0, 0)
    ) +
    labs(
      x = "Species",
      y = "Percentage contribution to significant tail-dependence",
      fill = NULL
    ) +
    theme_bw(base_size = 16) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black",
                                 size = label_size),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5,
                                face = "bold"),
      legend.text = element_text(color = "black")
    )
  
  return(p)
}

# Now visualization

# ============ 0-250km ===========

#nbin<-4
#target_dist_cat<-c(0,250)

gall_sig95CI32 <- visualize_spat_syn_Nyrs_threshold(yr_threshold = 32, 
                                                            target_dist_cat = c(0, 250), 
                                                            nbin = 4, 
                                                            siglevel = "95%CI",
                                                            label_size = 10)


gall_sig95CI36 <- visualize_spat_syn_Nyrs_threshold(yr_threshold = 36, 
                                                        target_dist_cat = c(0, 250), 
                                                        nbin = 4, 
                                                        siglevel = "95%CI",
                                                        label_size = 12)

gall_sig95CI40 <- visualize_spat_syn_Nyrs_threshold(yr_threshold = 40, 
                                                        target_dist_cat = c(0, 250), 
                                                        nbin = 4, 
                                                        siglevel = "95%CI",
                                                        label_size = 12)

p40<-plot_species_categories(yr_threshold = 40)
p36<-plot_species_categories(yr_threshold = 36)
p32<-plot_species_categories(yr_threshold = 32)

pdf(here("RESULTS/visualize_spat_syn_for_abund_0_250km_nbin_4_nogroup_sig95CI.pdf"), width = 12, height = 17) # Open a new pdf file
gridExtra::grid.arrange(p40, gall_sig95CI40,
                        layout_matrix = rbind(c(1),  # row 1 → p40
                                              c(2),  # row 2 → gall_sig95CI40
                                              c(2)   # row 3 → gall_sig95CI40 again (span)
                                              ))
dev.off() 

pdf(here("RESULTS/visualize_spat_syn_for_abund_0_250km_nbin_4_nogroup_sig95CI_32yr.pdf"), width = 12, height = 22) # Open a new pdf file
gridExtra::grid.arrange(p32, gall_sig95CI32,
                        layout_matrix = rbind(c(1),  
                                              c(2),  
                                              c(2)   
                        ))
dev.off() 

pdf(here("RESULTS/visualize_spat_syn_for_abund_0_250km_nbin_4_nogroup_sig95CI_36yr.pdf"), width = 12, height = 22) # Open a new pdf file
gridExtra::grid.arrange(p36, gall_sig95CI36,
                        layout_matrix = rbind(c(1),  
                                              c(2),  
                                              c(2)   
                        ))
dev.off() 






