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
boot32 <- readRDS("RESULTS/model_phylolm_sig95_0-250km_min32yr/model_tas5/model_est_phylolm_boot.RDS")
boot36 <- readRDS("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_tas5/model_est_phylolm_boot.RDS")
boot40 <- readRDS("RESULTS/model_phylolm_sig95_0-250km/model_tas5/model_est_phylolm_boot.RDS")

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


g1<-ggplot(data=dfsig.abs,
           aes(x=ftas5.sig,y=fab.sig,col=year), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+# existing overlays...
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.3)+
  xlab("Directional (Net) tail-dependent spatial synchrony,\n Temperature")+
  ylab("Directional (Net) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=year, fill=year), linetype = "dotted")+
  theme_bw(base_size = 18)+
  facet_wrap(~year)+
  stat_cor(aes(col=year,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(-1,-1,-1), label.y = c(0.8, 0.8,0.8), show.legend = F)+
  stat_regline_equation(label.x = c(-1,-1), label.y = c(0.9, 0.9,0.9),
                        aes(col=year),show.legend = F)+
#geom_abline(intercept = I, slope = s2, linetype = "dashed", data=tempo) +
  scale_color_manual(values=pal)+
  scale_fill_manual(values=pal)+
 theme(legend.position="none",panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(),legend.background=element_blank())

#g1


## ---- 1) Build prediction summary from ALL bootstrap draws (tidy pipeline) ----

x_grid <- seq(-1, 1, length.out = 200)

pred_df <-
  xall %>%
  # keep only the needed coefficients
  select(year, `(Intercept)`, ftas5.sig) %>%
  # give each bootstrap row an id (per year) so intercept/slope stay paired
  group_by(year) %>%
  mutate(draw = row_number()) %>%
  ungroup() %>%
  # make the x grid (same for all years)
  crossing(x = x_grid) %>%
  # predict for each draw across x
  mutate(y = `(Intercept)` + ftas5.sig * x) %>%
  # summarize across draws to median and 95% CI at each x
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
  x = c(-1, -1),       # x position
  y = c(1, 0.70),       # y positions
  label = c("Dotted + shade = OLS with 95%CI (no phylogeny)",
            "Solid + dashed = Phylogenetic bootstrap"),
  year = factor("40 years", levels = c("40 years","36 years","32 years"))
)
g1T<-g1 +
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
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))+
  geom_text(data = ann_text,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, size = 4)

#======================================================================================================
phylo_sum <- xall %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    draws = dplyr::n(),
    slope_med = median(ftas5.sig, na.rm = TRUE),
    slope_lo  = quantile(ftas5.sig, 0.025, na.rm = TRUE),
    slope_hi  = quantile(ftas5.sig, 0.975, na.rm = TRUE),
    p_boot    = 2 * pmin(mean(ftas5.sig >= 0), mean(ftas5.sig <= 0)), # two-sided
    .groups = "drop"
  )

phylo_sum
ann_p <- phylo_sum %>%
  dplyr::mutate(
    label = sprintf("Phylo slope 95%% CI: [%.2f, %.2f]\np_boot = %.3f",
                    slope_lo, slope_hi, p_boot),
    x = -1, y = 0.6
  )

g1_final <- g1T +                          # your current layered plot object
  geom_text(data = ann_p,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4)
g1_final

pdf(here("Results/plot_regression_0-250km_40to32years_netTaildep.pdf"),width=18, height=6)
print(g1_final)
dev.off()








