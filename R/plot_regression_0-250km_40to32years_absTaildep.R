#==================== plot across species ==============
rm(list=ls())
library(here)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(readxl)
nbin<-4

dfsig95.40<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km.csv",sep="")))# for 40 yrs
dfsig95.32<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km_min32yr.csv",sep="")))# for 32 yrs
dfsig95.36<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km_min36yr.csv",sep="")))# for 36 yrs


dfsig95.40.abs<-dfsig95.40%>%dplyr::select(abs.tot.td.pr5.sig,fpr5.sig,
                                           abs.tot.td.tas5.sig,ftas5.sig,
                                           abs.tot.td.ab.sig, fab.sig)%>%mutate(year="40")

dfsig95.32.abs<-dfsig95.32%>%dplyr::select(abs.tot.td.pr5.sig,fpr5.sig,
                                           abs.tot.td.tas5.sig,ftas5.sig,
                                           abs.tot.td.ab.sig, fab.sig)%>%mutate(year="32")

dfsig95.36.abs<-dfsig95.36%>%dplyr::select(abs.tot.td.pr5.sig,fpr5.sig,
                                           abs.tot.td.tas5.sig,ftas5.sig,
                                           abs.tot.td.ab.sig, fab.sig)%>%mutate(year="36")

dfsig.abs<-rbind(dfsig95.40.abs,dfsig95.36.abs,dfsig95.32.abs)
dfsig.abs$year<-as.factor(dfsig.abs$year)
dfsig.abs$year <- factor(dfsig.abs$year, 
                         levels = c("40", "36", "32"),
                         labels = c("40 years", "36 years", "32 years"))
#==============================================================================================================
boot32 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5_abs_td/model_est_phylolm_boot.RDS"))
boot36 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5_abs_td/model_est_phylolm_boot.RDS"))
boot40 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km/model_pr5_abs_td/model_est_phylolm_boot.RDS"))

x32<-boot32$bootstrap
x32<-as.data.frame(x32)
x32$year<-"32 years"

x36<-boot36$bootstrap
x36<-as.data.frame(x36)
x36$year<-"36 years"

x40<-boot40$bootstrap
x40<-as.data.frame(x40)
x40$year<-"40 years"

xall<-rbind(x40,x36,x32)
xall$year<-factor(xall$year, 
                  levels = c("40 years", "36 years", "32 years"),
                  labels = c("40 years", "36 years", "32 years"))

str(xall)
str(dfsig.abs)

## ---- 0) Ensure matching factor levels (robust) ----
levels_year <- c("40 years", "36 years", "32 years")
dfsig.abs$year <- factor(dfsig.abs$year, levels = levels_year)
xall$year      <- factor(xall$year,      levels = levels_year)

## Optional named palette (prevents color/level shuffles)
pal <- c("40 years" = "#000000", "36 years" = "#cc79a7", "32 years" = "#A6761D")

#============================================= for precipitation ===============================================
g1.wEq<-ggplot(data=dfsig.abs,
           aes(x=abs.tot.td.pr5.sig,y=abs.tot.td.ab.sig,col=year), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+# existing overlays...
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Precipitation")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=year, fill=year), linetype = "dotted")+
  theme_bw(base_size = 20)+
  facet_wrap(~year, scale="free_x")+
  stat_cor(aes(col=year,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(0,0,0), label.y = c(30,30,30), show.legend = F)+
  stat_regline_equation(label.x = c(0,0,0), label.y = c(31,31,31),
                        aes(col=year),show.legend = F)+
  scale_color_manual(values=pal)+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())

g1<-ggplot(data=dfsig.abs,
               aes(x=abs.tot.td.pr5.sig,y=abs.tot.td.ab.sig,col=year), 
               add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+# existing overlays...
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Precipitation")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=year, fill=year), linetype = "dotted")+
  theme_bw(base_size = 20)+
  facet_wrap(~year, scale="free")+
  scale_color_manual(values=pal)+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())


## ---- 1) Build prediction summary from ALL bootstrap draws (tidy pipeline) ----

# Get x range per year from the actual data
x_ranges <- dfsig.abs %>%
  group_by(year) %>%
  summarise(
    x_min = min(abs.tot.td.pr5.sig, na.rm = TRUE),
    x_max = max(abs.tot.td.pr5.sig, na.rm = TRUE),
    .groups = "drop"
  )

# Join x_ranges to bootstrap draws and generate x grids per year
pred_df <- xall %>%
  select(year, `(Intercept)`, abs.tot.td.pr5.sig) %>%
  group_by(year) %>%
  mutate(draw = row_number()) %>%
  ungroup() %>%
  inner_join(x_ranges, by = "year") %>%
  group_by(year) %>%
  # Create x grid per year based on actual observed x range
  do(crossing(., x = seq(.$x_min[1], .$x_max[1], length.out = 200))) %>%
  ungroup() %>%
  # Predict
  mutate(y = `(Intercept)` + abs.tot.td.pr5.sig * x) %>%
  group_by(year, x) %>%
  summarise(
    y50 = median(y, na.rm = TRUE),
    ylo = quantile(y, 0.025, na.rm = TRUE),
    yhi = quantile(y, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

## ---- 2) Overlay on your existing g1 (solid median + dotted CI; axes [-1, 1]) ----
# annotation data for the first facet only
ann_text <- data.frame(
  x = c(0, 0),       # x position
  y = c(20, 19),       # y positions
  label = c("Dotted + shade = OLS with 95%CI (no phylogeny)",
            "Solid + dashed = Phylogenetic bootstrap"),
  year = factor("40 years", levels = c("40 years","36 years","32 years"))
)
g1P.wEq<-g1.wEq +
  geom_line(data = pred_df,
            aes(x = x, y = y50, color = year),
            size = 1, inherit.aes = FALSE) +
  geom_line(data = pred_df,
            aes(x = x, y = ylo, color = year),
            linetype = "dashed", size=1, inherit.aes = FALSE) +
  geom_line(data = pred_df,
            aes(x = x, y = yhi, color = year),
            linetype = "dashed", size=1, inherit.aes = FALSE) +
  geom_text(data = ann_text,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, size = 4)

g1P<-g1 +
  geom_line(data = pred_df,
            aes(x = x, y = y50, color = year),
            size = 1, inherit.aes = FALSE) +
  geom_line(data = pred_df,
            aes(x = x, y = ylo, color = year),
            linetype = "dashed", size=1, inherit.aes = FALSE) +
  geom_line(data = pred_df,
            aes(x = x, y = yhi, color = year),
            linetype = "dashed", size=1, inherit.aes = FALSE)

tab32<-read.table(here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_pr5_abs_td/boot_summary.txt"), header=T, nrows=3)
tab36<-read.table(here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_pr5_abs_td/boot_summary.txt"), header=T, nrows=3)
tab40<-read.table(here("RESULTS/model_phylolm_sig95_0-250km/model_pr5_abs_td/boot_summary.txt"), header=T, nrows=3)

s40 <- tab40[tab40$term == "abs.tot.td.pr5.sig", ]
s36 <- tab36[tab36$term == "abs.tot.td.pr5.sig", ]
s32 <- tab32[tab32$term == "abs.tot.td.pr5.sig", ]

tab_all<-rbind(s40,s36,s32)
tab_all<-as.data.frame(tab_all)
tab_all$year<-factor(levels_year,levels = levels_year)

tab_all <- tab_all %>%
  mutate(
    label = sprintf("Slope Median = %.2f\nSlope 95%% CI: [%.2f, %.2f]\np_boot = %.3f",
                    slope_med, slope_lo, slope_hi, p_boot),
    x = 0,  # Set desired x position for annotation (adjust if needed)
    y = 25  # Set desired y position for annotation (adjust if needed)
)  

g1P.wEq <- g1P.wEq +
  geom_text(data = tab_all,
            aes(x = x, y = y, label = label, color = year),
            inherit.aes = FALSE,
            hjust = 0, size = 4, show.legend = FALSE)

#pdf(here("Results/plot_regression_0-250km_40to32years_absTaildep_P.pdf"),width=18, height=6)
print(g1P)
#dev.off()
#==========================================================================================================================
#================================================== now for temperature ==================================================
#==========================================================================================================================

boot32 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_tas5_abs_td/model_est_phylolm_boot.RDS"))
boot36 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_tas5_abs_td/model_est_phylolm_boot.RDS"))
boot40 <- readRDS(here("RESULTS/model_phylolm_sig95_0-250km/model_tas5_abs_td/model_est_phylolm_boot.RDS"))

x32<-boot32$bootstrap
x32<-as.data.frame(x32)
x32$year<-"32 years"

x36<-boot36$bootstrap
x36<-as.data.frame(x36)
x36$year<-"36 years"

x40<-boot40$bootstrap
x40<-as.data.frame(x40)
x40$year<-"40 years"

xall<-rbind(x40,x36,x32)
xall$year<-factor(xall$year, 
                  levels = c("40 years", "36 years", "32 years"),
                  labels = c("40 years", "36 years", "32 years"))

str(xall)
str(dfsig.abs)

## ---- 0) Ensure matching factor levels (robust) ----
levels_year <- c("40 years", "36 years", "32 years")
dfsig.abs$year <- factor(dfsig.abs$year, levels = levels_year)
xall$year      <- factor(xall$year,      levels = levels_year)

## Optional named palette (prevents color/level shuffles)
pal <- c("40 years" = "#000000", "36 years" = "#cc79a7", "32 years" = "#A6761D")

#============================================= for temperature ===============================================
g2.wEq<-ggplot(data=dfsig.abs,
           aes(x=abs.tot.td.tas5.sig,y=abs.tot.td.ab.sig,col=year), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+# existing overlays...
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Temperature")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=year, fill=year), linetype = "dotted")+
  theme_bw(base_size = 20)+
  facet_wrap(~year, scale="free_x")+
  stat_cor(aes(col=year,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(0,0,0), label.y = c(30,30,30), show.legend = F)+
  stat_regline_equation(label.x = c(0,0,0), label.y = c(31,31,31),
                        aes(col=year),show.legend = F)+
  scale_color_manual(values=pal)+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())

g2<-ggplot(data=dfsig.abs,
           aes(x=abs.tot.td.tas5.sig,y=abs.tot.td.ab.sig,col=year), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+# existing overlays...
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Temperature")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=year, fill=year), linetype = "dotted")+
  theme_bw(base_size = 20)+
  facet_wrap(~year, scale="free")+
  scale_color_manual(values=pal)+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())


## ---- 1) Build prediction summary from ALL bootstrap draws (tidy pipeline) ----

# Get x range per year from the actual data
x_ranges <- dfsig.abs %>%
  group_by(year) %>%
  summarise(
    x_min = min(abs.tot.td.tas5.sig, na.rm = TRUE),
    x_max = max(abs.tot.td.tas5.sig, na.rm = TRUE),
    .groups = "drop"
  )

# Join x_ranges to bootstrap draws and generate x grids per year
pred_df <- xall %>%
  select(year, `(Intercept)`, abs.tot.td.tas5.sig) %>%
  group_by(year) %>%
  mutate(draw = row_number()) %>%
  ungroup() %>%
  inner_join(x_ranges, by = "year") %>%
  group_by(year) %>%
  # Create x grid per year based on actual observed x range
  do(crossing(., x = seq(.$x_min[1], .$x_max[1], length.out = 200))) %>%
  ungroup() %>%
  # Predict
  mutate(y = `(Intercept)` + abs.tot.td.tas5.sig * x) %>%
  group_by(year, x) %>%
  summarise(
    y50 = median(y, na.rm = TRUE),
    ylo = quantile(y, 0.025, na.rm = TRUE),
    yhi = quantile(y, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

## ---- 2) Overlay on your existing g1 (solid median + dotted CI; axes [-1, 1]) ----
# annotation data for the first facet only
ann_text <- data.frame(
  x = c(0, 0),       # x position
  y = c(19, 18),       # y positions
  label = c("Dotted + shade = OLS with 95%CI (no phylogeny)",
            "Solid + dashed = Phylogenetic bootstrap"),
  year = factor("40 years", levels = c("40 years","36 years","32 years"))
)
g2T.wEq<-g2.wEq +
  geom_line(data = pred_df,
            aes(x = x, y = y50, color = year),
            size = 1, inherit.aes = FALSE) +
  geom_line(data = pred_df,
            aes(x = x, y = ylo, color = year),
            linetype = "dashed", size=1, inherit.aes = FALSE) +
  geom_line(data = pred_df,
            aes(x = x, y = yhi, color = year),
            linetype = "dashed", size=1, inherit.aes = FALSE) +
  scale_color_manual(values = pal) +
  geom_text(data = ann_text,
            aes(x = x, y = y, label = label),
           inherit.aes = FALSE, hjust = 0, size = 4)
  
g2T<-g2 +
    geom_line(data = pred_df,
              aes(x = x, y = y50, color = year),
              size = 1, inherit.aes = FALSE) +
    geom_line(data = pred_df,
              aes(x = x, y = ylo, color = year),
              linetype = "dashed", size=1, inherit.aes = FALSE) +
    geom_line(data = pred_df,
              aes(x = x, y = yhi, color = year),
              linetype = "dashed", size=1, inherit.aes = FALSE) +
  scale_color_manual(values = pal) 

tab32<-read.table(here("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_tas5_abs_td/boot_summary.txt"), header=T, nrows=3)
tab36<-read.table(here("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_tas5_abs_td/boot_summary.txt"), header=T, nrows=3)
tab40<-read.table(here("RESULTS/model_phylolm_sig95_0-250km/model_tas5_abs_td/boot_summary.txt"), header=T, nrows=3)

s40 <- tab40[tab40$term == "abs.tot.td.tas5.sig", ]
s36 <- tab36[tab36$term == "abs.tot.td.tas5.sig", ]
s32 <- tab32[tab32$term == "abs.tot.td.tas5.sig", ]

tab_all<-rbind(s40,s36,s32)
tab_all<-as.data.frame(tab_all)
tab_all$year<-factor(levels_year,levels = levels_year)

tab_all <- tab_all %>%
  mutate(
    label = sprintf("Slope Median = %.2f\nSlope 95%% CI: [%.2f, %.2f]\np_boot = %.3f",
                    slope_med, slope_lo, slope_hi, p_boot),
    x = 0,  # Set desired x position for annotation (adjust if needed)
    y = 25  # Set desired y position for annotation (adjust if needed)
  )  

g2T.wEq <- g2T.wEq +
  geom_text(data = tab_all,
            aes(x = x, y = y, label = label, color = year),
            inherit.aes = FALSE,
            hjust = 0, size = 4, show.legend = FALSE)

pdf(here("Results/plot_regression_0-250km_40to32years_absTaildep_PT_ref.pdf"),width=18, height=12)
grid.arrange(g1P.wEq,g2T.wEq, nrow=2)
dev.off()

pdf(here("Results/plot_regression_0-250km_40to32years_absTaildep_PT.pdf"),width=18, height=12)
grid.arrange(g1P,g2T, nrow=2)
dev.off()











