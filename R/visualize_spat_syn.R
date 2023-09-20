#rm(list=ls())
library(tidyverse)
library(gridExtra)
library(here)

visualize_spat_syn<-function(plotonly="LT",df1,df2){
  
  # select few variables
  df2<-df2%>%dplyr::select(AOU, ScientificName, totsites, ngoodsites,
                    SpecID, ORDER, Family, Genus, Species, 
                    PassNonPass, Diet.5Cat, Diet.Certainty, Strat.7Cat, ForStrat.SpecLevel,
                    nAbsentSitesAllyr, Nocturnal, IUCN_status)
  
  df<-left_join(df2,df1,by="AOU")# This is the dataframe we need to visualize
  
  df<-df%>%filter(nint!=nind)# exluding sp. where all interactions are indep.
  
  df$fL<-(df$nL/(df$nint))*100
  df$fU<-(df$nU/(df$nint))*100
  df$fLval<-(df$L/(df$L + abs(df$U)))*100
  df$fUval<-(abs(df$U)/(df$L + abs(df$U)))*100
  
  id<-which(is.na(df$fLval))
  df$fLval[id]<-0 # replacing NaN # NaN happens when all are either indep, or neg corr.
  df$fUval[id]<-0 # replacing NaN
  
  class(df$fLval)
  #plot(df$fL, df$fLval)
  
  table(df$Diet.5Cat)
  
  df$LandU<-df$fLval-df$fUval
  
  if(plotonly=="LT"){
    df<-df%>%filter(LandU>0)# only for positive tail dep. (lower tail) synchrony
  }
  
  if(plotonly=="UT"){
    df<-df%>%filter(LandU<0)# only for positive tail dep. (lower tail) synchrony
  }
  print(nrow(df))
  # arrange data in required format for stacked circular barplot
  dfsmall<-df%>%select(individual=AOU, group=Diet.5Cat, fLval, fUval)#%>%
  #              mutate(individual=paste(AOU,IUCN_status,sep=","))%>%select(-AOU,-IUCN_status)
  
  dfs2<-dfsmall%>%gather(key="observation", value="value", -c(1,2))
  
  #################################################################
  
  #> Loading required package: viridisLite
  data<-dfs2
  data$group<-substr(data$group,1,3)
  
  # Convert to factor
  data$group <- factor(data$group)
  
  # Transform data in a tidy format (long format)
  #data <- data %>% gather(key = "observation", value = "value", -c(1, 2))
  
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 2
  nObsType <- nlevels(as.factor(data$observation))
  to_add <- data.frame(matrix(NA, empty_bar * nlevels(data$group) * nObsType, ncol(data)))
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(data$group), each = empty_bar * nObsType)
  data <- rbind(data, to_add)
  data <- data %>% arrange(group, individual)
  data$id <- rep(seq(1, nrow(data) / nObsType), each = nObsType)
  
  # Get the name and the y position of each label
  label_data <- data %>%
    group_by(id, individual) %>%
    summarize(tot = sum(value))
  #> `summarise()` has grouped output by 'id'. You can override using the `.groups`
  #> argument.
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)
  
  # prepare a data frame for base lines
  base_data <- data %>%
    group_by(group) %>%
    summarize(start = min(id), end = max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title = mean(c(start, end)))
  
  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1, ]
  
  # Make the plot
  g1<-ggplot(data) +
    
    # Add the stacked bar
    geom_bar(aes(x = as.factor(id), y = value, fill = observation), stat = "identity", alpha = 0.8) +
    #scale_fill_viridis(discrete = TRUE) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha = 1, size = 0.3, inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "black", alpha = 1, size = 0.3, inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "black", alpha = 1, size = 0.3, inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "black", alpha = 1, size = 0.3, inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "black", alpha = 1, size = 0.3, inherit.aes = FALSE) +
    
    # Add text showing the value of each 100/75/50/25 lines
    ggplot2::annotate("text", x = rep(max(data$id), 5), y = c(0, 25, 50, 75, 100), 
                      label = c("0", "25", "50", "75", "100"), 
                      color = "black", size = 3, angle = 0, fontface = "bold", hjust = 0.8) +
    ylim(-150, max(label_data$tot, na.rm = T)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm")
    ) +
    coord_polar() +
    
    # Add labels on top of each bar
    geom_text(data = label_data, aes(x = id, y = tot -0.5, 
                                     label = individual, 
                                     hjust = hjust), color = "gray", 
              fontface = "bold", alpha = 0.6, size = 4, 
              angle = label_data$angle, inherit.aes = FALSE) +
    
    # Add base line information
    geom_segment(data = base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha = 0.8, size = 1, inherit.aes = FALSE) +
    geom_text(data = base_data, aes(x = title, y = -10, label = group), hjust = c(1, 1.0, 0, 0,0.8), colour = "black", alpha = 0.8, size = 3, fontface = "bold", angle=30,inherit.aes = FALSE)
  #> Warning: Removed 24 rows containing missing values (position_stack).
  #> Warning: Removed 9 rows containing missing values (geom_text).
  return(g1)
}

# Now visualization

# ============ 0-250km ===========
plotonly<-"LT"
df1<-read.csv(here("RESULTS/summary_spat_syn_for_abund_0_250km.csv"))
df2<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
gLT<-visualize_spat_syn(df1=df1,df2=df2,plotonly=plotonly)

plotonly<-"UT"
df1<-read.csv(here("RESULTS/summary_spat_syn_for_abund_0_250km.csv"))
df2<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
gUT<-visualize_spat_syn(df1=df1,df2=df2,plotonly=plotonly)

# later in your plot you could add info about IUCN status in inkscape

pdf(here("RESULTS/visualize_spat_syn_for_abund_0_250km.pdf"), width = 14, height = 7) # Open a new pdf file
grid.arrange(gLT, gUT, nrow=1) # Write the grid.arrange in the file
dev.off() 


#df263<-rbind(dfl,dfu)
#df263<-df263[,c(1,2)]
#write.csv(df263,here("RESULTS/species263_127LT_136UT.csv"),row.names = F)
#df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_with_optimal_biovar.csv"))
#id254<-which(df263$AOU%in%df$AOU)
#df263$AOU[setdiff(1:263,id254)] # not included in phylotree

