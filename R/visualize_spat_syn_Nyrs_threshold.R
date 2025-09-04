# --- Packages ---
rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)

#--------------------

visualize_spat_syn_Nyrs_threshold<-function(yr_threshold=40,target_dist_cat=c(0,250), nbin=4, siglevel="95%CI"){
  
  if(yr_threshold==40){
    df1<-readRDS(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,
                            "_corlmcoru_sigres_summary_",target_dist_cat[1],"-",
                            target_dist_cat[2],"Km.RDS",sep="")))
  }else{
    df1<-readRDS(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,
                            "_corlmcoru_sigres_summary_",target_dist_cat[1],"-",
                            target_dist_cat[2],"Km_min",yr_threshold,"yr.RDS",sep="")))
  }
  
  if(siglevel=="75%CI"){
    df1<-df1%>%filter(Lsig75ab!=0 | Usig75ab!=0)
    df1<-df1%>%dplyr::select(AOU,Lsigab=Lsig75ab, Usigab=Usig75ab)
  }else{
    df1<-df1%>%filter(Lsig95ab!=0 | Usig95ab!=0)
    df1<-df1%>%dplyr::select(AOU,Lsigab=Lsig95ab, Usigab=Usig95ab)
  }
  
  
  df1$fLval<-(df1$Lsigab/(df1$Lsigab + abs(df1$Usigab)))*100
  df1$fUval<-(abs(df1$Usigab)/(df1$Lsigab + abs(df1$Usigab)))*100
  
  # arrange data in required format for stacked circular barplot
  dfsmall<-df1%>%dplyr::select(individual=AOU, fLval, fUval)
  data<-dfsmall%>%gather(key="observation", value="value", -c(1))
  #################################################################
  
  
  # --- TUNABLES ---
  inner_radius <- 120       # bigger = bigger white hole
  order_by_total <- TRUE    # TRUE: order bars by total height; FALSE: order by individual
  label_size <- 3
  
  # --- PREP DATA (robust) ---
  df <- data %>%
    mutate(
      individual = as.character(individual),        # ensure consistent type
      observation = as.character(observation)
    ) %>%
    filter(observation %in% c("fLval","fUval")) %>%
    mutate(observation = factor(observation, levels = c("fLval","fUval"))) %>%
    complete(individual, observation, fill = list(value = 0)) %>%
    group_by(individual, observation) %>%
    summarize(value = mean(value, na.rm = TRUE), .groups = "drop")  # use sum() if preferred
  
  # totals and desired order
  totals <- df %>%
    group_by(individual) %>%
    summarize(total = sum(value, na.rm = TRUE), .groups = "drop")
  
  ord_vec <- if (order_by_total) {
    totals %>% arrange(desc(total)) %>% pull(individual)
  } else {
    totals %>% arrange(individual) %>% pull(individual)
  }
  
  # add an explicit bar_id via safe matching
  dfp <- df %>%
    mutate(bar_id = match(individual, ord_vec)) %>%
    arrange(bar_id, observation)
  
  # label data (one row per bar, at top of stack)
  label_data <- totals %>%
    mutate(bar_id = match(individual, ord_vec)) %>%
    arrange(bar_id)
  
  nbar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$bar_id - 0.5) / nbar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)
  
  ymax <- max(label_data$total, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1  # safety
  
  # --- PLOT: circular stacked donut ---
  ref_vals <- c(0, 25, 50, 75, 100)
  
  p <- ggplot(dfp, aes(x = as.factor(bar_id), y = value, fill = observation)) +
    # --- background reference rings ---
    geom_hline(yintercept = ref_vals, colour = "grey60", size = 0.5) +
    
    # --- stacked bars on top ---
    geom_col(width = 0.9, alpha = 0.95) +
    coord_polar() +
    ylim(-inner_radius, ymax * 1.15) +
    
    # outer labels for individuals
    geom_text(
      data = label_data,
      aes(x = bar_id, y = total + 0.03 * ymax,
          label = individual, hjust = hjust),
      inherit.aes = FALSE,
      angle = label_data$angle,
      size = label_size,
      fontface = "bold"
    ) +
    
    # --- reference labels drawn AFTER lines ---
    annotate("text",
             x = nrow(label_data) + 3,   # position just outside bars
             y = ref_vals,
             label = ref_vals,
             colour = "grey60",
             size = 3.5,
             fontface = "bold",
             hjust = 1,
             vjust = -0.3) +               # move slightly above the line
    
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-0.5, 4), "cm")
    )
  
  return(p)
  
}

# Now visualization

# ============ 0-250km ===========

#nbin<-4
#target_dist_cat<-c(0,250)

gall_sig75CI<-visualize_spat_syn_Nyrs_threshold(yr_threshold = 32, target_dist_cat= c(0,250), nbin=4, siglevel="75%CI")
gall_sig95CI<-visualize_spat_syn_Nyrs_threshold(yr_threshold = 32, target_dist_cat= c(0,250), nbin=4, siglevel="95%CI")

pdf(here("RESULTS/visualize_spat_syn_for_abund_0_250km_nbin_4_nogroup_sig75_and_95CI_min32yr.pdf"), width = 16, height = 10) # Open a new pdf file
gridExtra::grid.arrange(gall_sig75CI,gall_sig95CI, ncol=2)# Write the grid.arrange in the file
dev.off() 

gall_sig75CI<-visualize_spat_syn_Nyrs_threshold(yr_threshold = 36, target_dist_cat= c(0,250), nbin=4, siglevel="75%CI")
gall_sig95CI<-visualize_spat_syn_Nyrs_threshold(yr_threshold = 36, target_dist_cat= c(0,250), nbin=4, siglevel="95%CI")

pdf(here("RESULTS/visualize_spat_syn_for_abund_0_250km_nbin_4_nogroup_sig75_and_95CI_min36yr.pdf"), width = 16, height = 10) # Open a new pdf file
gridExtra::grid.arrange(gall_sig75CI,gall_sig95CI, ncol=2)# Write the grid.arrange in the file
dev.off() 

gall_sig75CI<-visualize_spat_syn_Nyrs_threshold(yr_threshold = 40, target_dist_cat= c(0,250), nbin=4, siglevel="75%CI")
gall_sig95CI<-visualize_spat_syn_Nyrs_threshold(yr_threshold = 40, target_dist_cat= c(0,250), nbin=4, siglevel="95%CI")

pdf(here("RESULTS/visualize_spat_syn_for_abund_0_250km_nbin_4_nogroup_sig75_and_95CI.pdf"), width = 16, height = 10) # Open a new pdf file
gridExtra::grid.arrange(gall_sig75CI,gall_sig95CI, ncol=2)# Write the grid.arrange in the file
dev.off() 







