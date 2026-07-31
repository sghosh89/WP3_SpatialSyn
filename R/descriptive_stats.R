rm(list=ls())
library(here)
library(tidyverse)
library(purrr)
library(PupillometryR)
library(patchwork)

data_40<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))
data_32<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min32yr.csv"))
data_36<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min36yr.csv"))

abund_array40<-readRDS(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_array.RDS"))
dim(abund_array40)# 41 by 161 by 498
dimnames(abund_array40)[[2]]

abund_array36<-readRDS(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_array_min36yr.RDS"))
dim(abund_array36) #41 by 575 by 564

abund_array32<-readRDS(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_array_min32yr.RDS"))
dim(abund_array32)#41 by 1018 by 615

# species dimension names of array must be AOU codes
# e.g. dimnames(abund_array40)[[3]]


#===============================
# 1. Combine datasets
#===============================

occ_df <- bind_rows(
  data_40 %>%
    mutate(Threshold = "40 years"),
  
  data_36 %>%
    mutate(Threshold = "36 years"),
  
  data_32 %>%
    mutate(Threshold = "32 years")
) %>%
  mutate(
    Threshold = factor(
      Threshold,
      levels = c("40 years", "36 years", "32 years")
    ),
    Occupancy = 100 * ngoodsites / totsites
  )

#===============================
# 2. Table S1
#===============================

table_s1 <- occ_df %>%
  group_by(Threshold) %>%
  summarise(
    `Total sites` = unique(totsites),
    `Species retained` = n(),
    
    `Median no. sites/species` = median(ngoodsites),
    `IQR no. sites/species` = paste0(
      quantile(ngoodsites, 0.25, type=1, na.rm=TRUE), "–",
      quantile(ngoodsites, 0.75, type=1, na.rm=TRUE)
    ),
    `Range no. sites/species` = paste0(
      min(ngoodsites), "–", max(ngoodsites)
    ),
    
    `Median occupancy (%)` = round(median(Occupancy), 2),
    `IQR occupancy (%)` = paste0(
      round(quantile(Occupancy, 0.25), 2), "–",
      round(quantile(Occupancy, 0.75), 2)
    ),
    `Range occupancy (%)` = paste0(
      round(min(Occupancy), 2), "–",
      round(max(Occupancy), 2)
    ),
    .groups = "drop"
  )

table_s1

#"Table S1. Species site occupancy across minimum time-series length thresholds."


write.csv(
  table_s1,
  "Results/Table_species_occupancy_summary.csv",
  row.names = FALSE
)

#===============================
# 3. Horizontal raincloud plot
#===============================

p_occ <- ggplot(
  occ_df,
  aes(x = Threshold, y = Occupancy, fill = Threshold, colour = Threshold)
) +
  
  geom_flat_violin(
    position = position_nudge(x = 0.18),
    trim = FALSE,
    alpha = 0.45,
    colour = NA
  ) +
  
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 1.2,
    alpha = 0.45
  ) +
  
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = NA,
    colour = "black"
  ) +
  
  coord_flip() +
  
  labs(
    x = "Minimum time-series length",
    y = "Site-coverage (%)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black")
  )

p_occ

ggsave(
  "Results/species_occupancy_raincloud.pdf",
  p_occ,
  width = 7,
  height = 6,
  dpi = 600
)

#===============================
xroutes<-read.csv(here("DATA/for_BBS/wrangled_data/uRID_lonlat_stratum.csv"))
#df_comsp<-unique(data_32$AOU,data_36$AOU,data_40$AOU)
#df<-data.frame(AOU=df_comsp)


#=========================================================
# Function
#=========================================================

summarize_climate <- function(prdat,
                              tempdat,
                              abund_array,
                              threshold_label){
  
  #all(colnames(prdat)==colnames(tempdat))==TRUE
  
  #---------------
  sitenames<-dimnames(abund_array)[[2]]
  entiresites<-colnames(prdat)
    
  id<-which(entiresites%in%sitenames)
  prdat_sel<-prdat[,id]
  tempdat_sel<-tempdat[,id]
  
  prdat_sel<-as.data.frame(prdat_sel)
  tempdat_sel<-as.data.frame(tempdat_sel)
  #-----------------------------
  # Long format
  #-----------------------------
  pr_long <- prdat_sel %>%
    mutate(Year = rownames(prdat_sel)) %>%
    pivot_longer(
      cols = -Year,
      names_to = "Site",
      values_to = "Precipitation"
    )
  
  temp_long <- tempdat_sel %>%
    mutate(Year = rownames(tempdat_sel)) %>%
    pivot_longer(
      cols = -Year,
      names_to = "Site",
      values_to = "Temperature"
    )
  
  clim_long <- left_join(
    temp_long,
    pr_long,
    by = c("Year", "Site")
  ) %>%
    mutate(
      Year = as.numeric(Year),
      Threshold = threshold_label,
      
      # Unit conversion
      Temperature_C = (Temperature/10) -273.15,
      Precipitation_mm = Precipitation/10
    )
  
  #-----------------------------
  # Summary table
  #-----------------------------
  summary_table <- clim_long %>%
    summarise(
      Threshold = threshold_label,
      Sites = n_distinct(Site),
      Years = n_distinct(Year),
      SiteYears = n(),
      
      TempMedian = median(Temperature_C, na.rm = TRUE),
      TempQ1 = quantile(Temperature_C, 0.25, na.rm = TRUE),
      TempQ3 = quantile(Temperature_C, 0.75, na.rm = TRUE),
      TempMin = min(Temperature_C, na.rm = TRUE),
      TempMax = max(Temperature_C, na.rm = TRUE),
      
      PrMedian = median(Precipitation_mm, na.rm = TRUE),
      PrQ1 = quantile(Precipitation_mm, 0.25, na.rm = TRUE),
      PrQ3 = quantile(Precipitation_mm, 0.75, na.rm = TRUE),
      PrMin = min(Precipitation_mm, na.rm = TRUE),
      PrMax = max(Precipitation_mm, na.rm = TRUE)
    )
  
  #-----------------------------
  # Annual means
  #-----------------------------
  yearly_summary <- clim_long %>%
    group_by(Year) %>%
    summarise(
      mean_temp = mean(Temperature_C, na.rm = TRUE),
      sd_temp = sd(Temperature_C, na.rm = TRUE),
      mean_pr = mean(Precipitation_mm, na.rm = TRUE),
      sd_pr = sd(Precipitation_mm, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(
    long = clim_long,
    summary = summary_table,
    yearly = yearly_summary
  ))
}

prdat<-readRDS(here("RESULTS/year_by_site_pr_avgAprtoAug.RDS"))
tempdat<-readRDS(here("RESULTS/year_by_site_tas_avgAprtoAug.RDS"))

clim40 <- summarize_climate(prdat, tempdat, abund_array=abund_array40, threshold_label="40 years")

clim36 <- summarize_climate(prdat, tempdat, abund_array=abund_array36, threshold_label="36 years")

clim32 <- summarize_climate(prdat, tempdat, abund_array=abund_array32, threshold_label="32 years")

climate_summary_table <- bind_rows(clim40$summary, clim36$summary, clim32$summary)

write.csv(
  climate_summary_table,
  "Results/Table_climate_summary.csv",
  row.names = FALSE
)

temp_ylim <- range(
  clim40$yearly$mean_temp + c(-clim40$yearly$sd_temp, clim40$yearly$sd_temp),
  clim36$yearly$mean_temp + c(-clim36$yearly$sd_temp, clim36$yearly$sd_temp),
  clim32$yearly$mean_temp + c(-clim32$yearly$sd_temp, clim32$yearly$sd_temp),
  na.rm = TRUE
)

pr_ylim <- range(
  clim40$yearly$mean_pr + c(-clim40$yearly$sd_pr, clim40$yearly$sd_pr),
  clim36$yearly$mean_pr + c(-clim36$yearly$sd_pr, clim36$yearly$sd_pr),
  clim32$yearly$mean_pr + c(-clim32$yearly$sd_pr, clim32$yearly$sd_pr),
  na.rm = TRUE
)
#-------------------------------------------------------
# Function to make annual climate plots
#-------------------------------------------------------

plot_climate <- function(df, threshold){
  
  p_temp <-
    ggplot(df, aes(Year, mean_temp)) +
    geom_ribbon(aes(ymin = mean_temp-sd_temp,
                    ymax = mean_temp+sd_temp),
                fill = "firebrick",
                alpha = 0.20) +
    geom_line(colour = "firebrick", linewidth = 0.7) +
    scale_x_continuous(
      breaks = c(1979, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2019),
      limits = c(1979, 2019),
      expand = c(0, 0)
    )+coord_cartesian(ylim = temp_ylim)+
    labs(
      title = threshold,
      x = NULL,
      y = expression(Temperature~(degree*C))
    ) +
    theme_classic(base_size = 13)
  
  p_pr <-
    ggplot(df, aes(Year, mean_pr)) +
    geom_ribbon(aes(ymin = mean_pr-sd_pr,
                    ymax = mean_pr+sd_pr),
                fill = "dodgerblue3",
                alpha = 0.20) +
    geom_line(colour = "dodgerblue4", linewidth = 0.7) +
    scale_x_continuous(
      breaks = c(1979, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2019),
      limits = c(1979, 2019),
      expand = c(0, 0)
    )+coord_cartesian(ylim = pr_ylim)+
    labs(
      x = "Year",
      y = "Precipitation (mm)"
    ) +
    theme_classic(base_size = 13)
  
  list(temp = p_temp,
       pr = p_pr)
  
}
p40 <- plot_climate(clim40$yearly,"minimum 40 years sampled")
p36 <- plot_climate(clim36$yearly,"minimum 36 years sampled")
p32 <- plot_climate(clim32$yearly,"minimum 32 years sampled")

climate_fig <-
  
  (p40$temp | p36$temp | p32$temp) /
  (p40$pr   | p36$pr   | p32$pr)

climate_fig

ggsave(
  "Results/Fig_climate_timeseries.pdf",
  climate_fig,
  width = 16,
  height = 6,
  dpi = 600
)
