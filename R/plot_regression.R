# run later below code
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

# plot for absolute tail-dep synchrony
gab1<-ggplot(data=dfsig.abs,
           aes(x=abs.tot.td.pr5.sig,y=abs.tot.td.ab.sig,col=year), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.3)+
  #ylim(0,5)+
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Precipitation")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=year, fill=year))+
  theme_bw()+
  stat_cor(aes(col=year,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(0.05,0.05,0.5), label.y = c(30, 30, 30), show.legend = F)+
  stat_regline_equation(label.x = c(0.05,0.05, 0.5), label.y = c(32,32, 32),
                        aes(col=year),show.legend = F)+
  facet_wrap(~year, scales = "free_x")+
  scale_color_manual(values=c("#000000", "#cc79a7","#A6761D"))+
  scale_fill_manual(values=c("#000000", "#cc79a7","#A6761D"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())
gab1

gab2<-ggplot(data=dfsig.abs,
             aes(x=abs.tot.td.tas5.sig,y=abs.tot.td.ab.sig,col=year), 
             add = "reg.line")+
  geom_point(pch=19, alpha=0.3)+
  #ylim(0,5)+
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Temperature")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=year, fill=year))+
  theme_bw()+
  stat_cor(aes(col=year,
             label = paste(..r.label..,..rr.label.., ..p.label.., 
                           sep = "*`,`~")),
         label.x = c(0.05,0.05,0.5), label.y = c(30, 30,30), show.legend = F)+
  stat_regline_equation(label.x = c(0.05,0.05,0.5), label.y = c(32,32,32),
                        aes(col=year),show.legend = F)+
  facet_wrap(~year, scales = "free_x")+
  scale_color_manual(values=c("#000000", "#cc79a7","#A6761D"))+
  scale_fill_manual(values=c("#000000", "#cc79a7","#A6761D"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())
gab2

pdf(here("Results/plot_regression_0-250km_for_abs_taildep_40to32years.pdf"),width=12.5, height=10)
grid.arrange(gab1,gab2,nrow=2)
dev.off()

#----------------------
slope_tab_dir<-read_excel("RESULTS/table_summary_phylolm.xlsx")

#g1<-ggplot(data=dfsig.abs,
#           aes(x=fpr5.sig,y=fab.sig,col=year), 
#           add = "reg.line")+
#  geom_point(pch=19, alpha=0.5)+
#  ylim(-1,1)+
#  xlab("Directional (Net) tail-dependent spatial synchrony,\n Precipitation")+
#  ylab("Directional (Net) tail-dependent spatial synchrony,\n Abundance")+
#  geom_smooth(method="lm", se=T,aes(col=year, fill=year),linetype = "dashed")+
#  theme_bw()+
#  facet_wrap(~year)+
#  stat_cor(aes(col=year,
#               label = paste(..r.label..,..rr.label.., ..p.label.., 
#                             sep = "*`,`~")),
#           label.x = c(-1,-1,-1), label.y = c(0.8, 0.8,0.8), show.legend = F)+
#  stat_regline_equation(label.x = c(-1,-1), label.y = c(0.9, 0.9,0.9),
#                        aes(col=year),show.legend = F)+
#  scale_color_manual(values=c("#000000", "#cc79a7","#A6761D"))+
#  scale_fill_manual(values=c("#000000", "#cc79a7","#A6761D"))+
#  theme(legend.position="none",panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),legend.background=element_blank())
#g1

tempo<-slope_tab_dir%>%filter(`Model coefficients`%in%c("I","s1") & model==4)# net tail-dep synchrony in temperature
to_num <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%                  # trim ends
    stringr::str_replace_all("[[:space:]]+", "") %>%  # remove internal spaces
    stringr::str_replace_all("\u2212", "-") %>%       # unicode minus → ascii minus
    readr::parse_number()                    # keeps decimals, scientific notation
}
tempo <- tempo %>%
  mutate(coef=`Model coefficients`,
    est = to_num(`Estimate`),
    lo = to_num(`lowerboot95%CI`),
    hi = to_num(`upperboot95%CI`),
    pvalue = to_num(`p-value`),
    year=as.factor(year)
  )
str(tempo)
tempo<-tempo%>%select(coef,est,lo,hi,pvalue,year)
# --- clean & reshape coef# --- clean & reshape `tempo` to wide (one row per year) ---
tempo_wide <- tempo %>%
  select(year, coef, est, lo, hi, pvalue) %>%
  pivot_wider(
    id_cols = year,
    names_from = coef,
    values_from = c(est, lo, hi, pvalue),
    names_glue = "{.value}_{coef}"
  ) 
# --- keep CI cols and s1 p-value for lines/labels ---
tempo_lines <- tempo_wide %>%
  mutate(
    year_num = as.numeric(as.character(year)),
    year = factor(paste0(year_num, " years"),
                  levels = levels(dfsig.abs$year))
  ) %>%
  select(year, est_I, est_s1, lo_I, lo_s1, hi_I, hi_s1, pvalue_s1)

# --- per-year observed x-range, clamped to [-1, 1] ---
xr <- dfsig.abs %>%
  group_by(year) %>%
  summarise(xmin = min(ftas5.sig, na.rm = TRUE),
            xmax = max(ftas5.sig, na.rm = TRUE), .groups = "drop") %>%
  mutate(xmin = pmax(xmin, -1), xmax = pmin(xmax,  1))

x_grid <- xr %>%
  rowwise() %>%
  mutate(x = list(if (xmin < xmax) seq(xmin, xmax, length.out = 200) else xmin)) %>%
  tidyr::unnest(x)

# --- central line and CI boundaries (conservative from coeff bounds) ---
# ensure bounds are ordered (just in case)
tempo_lines <- tempo_lines %>%
  mutate(
    lo_I  = pmin(lo_I, hi_I),  hi_I  = pmax(lo_I, hi_I),
    lo_s1 = pmin(lo_s1, hi_s1), hi_s1 = pmax(lo_s1, hi_s1)
  )

pred_df <- x_grid %>%
  left_join(tempo_lines, by = "year") %>%
  mutate(
    y_est = est_I + est_s1 * x,
    # corner envelope (guarantees y_lo <= y_est <= y_hi for all x)
    c1 = lo_I + lo_s1 * x,
    c2 = lo_I + hi_s1 * x,
    c3 = hi_I + lo_s1 * x,
    c4 = hi_I + hi_s1 * x,
    y_lo = pmin(c1, c2, c3, c4),
    y_hi = pmax(c1, c2, c3, c4)
)

# --- s1 p-value label (only) ---
p_labels <- tempo_lines %>%
  transmute(
    year,
    x = -0.95, y = 0.90,
    label = paste0("slope[s1]~p==", formatC(pvalue_s1, format = "g", digits = 3))
  )

# --- your plot: remove geom_abline(...) and use these three geom_line() calls ---
g1 <- ggplot(dfsig.abs, aes(x = ftas5.sig, y = fab.sig, col = year)) +
  geom_point(pch = 19, alpha = 0.5) +
  xlab("Directional (Net) tail-dependent spatial synchrony,\n Temperature") +
  ylab("Directional (Net) tail-dependent spatial synchrony,\n Abundance") +
  geom_smooth(method = "lm", se = TRUE, aes(fill = year), linetype = "dashed") +
  facet_wrap(~year) +
  stat_cor(aes(col = year,
               label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = c(-1,-1,-1), label.y = c(0.8, 0.8, 0.8), show.legend = FALSE) +
  stat_regline_equation(label.x = c(-1,-1), label.y = c(0.9, 0.9, 0.9),
                        aes(col = year), show.legend = FALSE) +
  # --- central solid line
  geom_line(data=pred_df, aes(x=x, y=y_est, color=year),
            inherit.aes=FALSE, linetype="solid",  linewidth=1) +
  geom_line(data=pred_df, aes(x=x, y=y_lo,  color=year),
            inherit.aes=FALSE, linetype="dotted", linewidth=0.7) +
  geom_line(data=pred_df, aes(x=x, y=y_hi,  color=year),
            inherit.aes=FALSE, linetype="dotted", linewidth=0.7) +
  # s1 p-value only
  geom_text(
    data = p_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, parse = TRUE, hjust = 0, vjust = 1, size = 3.5
  ) +
  scale_color_manual(values = c("#000000", "#cc79a7", "#A6761D")) +
  scale_fill_manual(values  = c("#000000", "#cc79a7", "#A6761D")) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank())

g1



g2<-ggplot(data=dfsig.abs,
           aes(x=ftas5.sig,y=fab.sig,col=year), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+
  xlab("Directional (Net) tail-dependent spatial synchrony,\n Temperature")+
  ylab("Directional (Net) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=year, fill=year), linetype = "dashed")+
  theme_bw()+
  facet_wrap(~year)+
  stat_cor(aes(col=year,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(-1,-1,-1), label.y = c(0.8, 0.8,0.8), show.legend = F)+
  stat_regline_equation(label.x = c(-1,-1), label.y = c(0.9, 0.9,0.9),
                        aes(col=year),show.legend = F)+
  #geom_abline(intercept = I, slope = s2, linetype = "dashed", data=tempo) +
  scale_color_manual(values=c("#000000", "#cc79a7","#A6761D"))+
  scale_fill_manual(values=c("#000000", "#cc79a7","#A6761D"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())


g2+
  geom_ribbon(
    data = pred_df,
    aes(x = x, ymin = ymin, ymax = ymax, fill = year),
    inherit.aes = FALSE, alpha = 0.15
  ) +
  geom_line(
    data = pred_df,
    aes(x = x, y = y_est, color = year),
    inherit.aes = FALSE, linewidth = 0.7, linetype = "dashed"
  )



pdf(here("Results/plot_regression_0-250km_40to32years.pdf"),width=12.5, height=10)
grid.arrange(g1,g2,nrow=2)
dev.off()




library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(rlang)

# -------- 1) Load bootstrap results (only 36 years) --------
boot <- readRDS("RESULTS/model_phylolm_sig95_0-250km_min36yr/model_tas5/model_est_phylolm_boot.RDS")

# Turn whatever came out into a data frame of coefficients per bootstrap draw
# Handles common cases: list of fitted models OR a data frame/matrix of coefs
if (is.list(boot) && !inherits(boot, "data.frame")) {
  # list of model objects
  coef_mat <- sapply(boot, function(m) coef(m))
  boot_df  <- as_tibble(t(coef_mat))
} else {
  boot_df <- as_tibble(boot)
}

# -------- 2) Pick the intercept and the coefficient for your x variable --------
xvar <- "ftas5.sig"  # <- your x variable name in dfsig.abs

cn <- make.names(names(boot_df))        # normalize names in case of spaces/parens
names(boot_df) <- cn

# find intercept column
int_name <- cn[cn %in% c(".Intercept.","Intercept","I","(Intercept)")]
if (length(int_name) == 0) int_name <- cn[grepl("Intercept|^I$", cn, ignore.case = TRUE)][1]

# find coefficient for xvar
x_name  <- cn[grepl(paste0("\\b", make.names(xvar), "\\b"), cn)]
if (length(x_name) == 0) x_name <- cn[grepl("s1|beta1|b1", cn, ignore.case = TRUE)][1]

stopifnot(length(int_name) == 1, length(x_name) == 1)

draws <- boot_df %>%
  transmute(I  = .data[[int_name]],
            s1 = .data[[x_name]]) %>%
  filter(is.finite(I), is.finite(s1))

# Optional: if you wanted to hold HWI at its mean rather than ignore it:
# hwi_name <- cn[grepl("HWI|s2|beta2|b2", cn, ignore.case = TRUE)][1]
# if (!is.na(hwi_name)) {
#   HWI_mean_36 <- dfsig.abs %>% filter(year == "36 years") %>%
#                   summarize(m = mean(HWI, na.rm = TRUE)) %>% pull(m)
#   draws <- boot_df %>%
#     transmute(I  = .data[[int_name]] + .data[[hwi_name]] * HWI_mean_36,
#               s1 = .data[[x_name]]) %>%
#     filter(is.finite(I), is.finite(s1))
# }

# -------- 3) Build x grid from observed 36-years data, clamped to [-1,1] --------
df36 <- dfsig.abs %>% filter(year == "36 years")
x_min <- max(-1, min(df36$ftas5.sig, na.rm = TRUE))
x_max <- min( 1, max(df36$ftas5.sig, na.rm = TRUE))
x_seq <- seq(x_min, x_max, length.out = 200)

# pointwise bootstrap quantiles for I + s1 * x
y_q <- sapply(x_seq, function(xx) {
  y <- draws$I + draws$s1 * xx
  quantile(y, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
})

pred_boot_36 <- tibble(
  year  = factor("36 years", levels = levels(dfsig.abs$year)),
  x     = x_seq,
  y_lo  = y_q[1,],
  y_med = y_q[2,],
  y_hi  = y_q[3,]
)

# -------- 4) Add to your existing plot (or build anew) --------
# If you already have g1, just add these three layers.
# If you were drawing the earlier "tempo" envelope, you may want to remove it now.

g1 +
  geom_line(
    data = pred_boot_36,
    aes(x = x, y = y_med, color = year),
    inherit.aes = FALSE, linewidth = 1
  ) +
  geom_line(
    data = pred_boot_36,
    aes(x = x, y = y_lo, color = year),
    inherit.aes = FALSE, linetype = "dotted", linewidth = 0.7
  ) +
  geom_line(
    data = pred_boot_36,
    aes(x = x, y = y_hi, color = year),
    inherit.aes = FALSE, linetype = "dotted", linewidth = 0.7
  ) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))
