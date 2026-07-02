############################################################
# Conceptual figure: copula plot and corresponding time series
# Tail-dependent points only are colored in B1/C1 and B2/C2
############################################################

library(VineCopula)
library(copula)
library(here)

source(here("R/ExtremeTailDep.R"))

set.seed(123)

#===========================================================
# Simulate data
#===========================================================

np <- 1000000
N  <- 40
time <- 1:N

# Copula with extreme right tail
xtc_r <- retd(n = np, d = 2, rl = 1, mn = 0, sdev = 1)
xtc_r <- pnorm(xtc_r)

# Convert to lower-tail dependent version
xtc_l <- 1 - xtc_r

# Frank copula: no tail dependence
cf <- iRho(frankCopula(), 0.875)
f_cop <- BiCopSim(N, family = 5, par = cf)

#===========================================================
# Define tail thresholds
#===========================================================

q_tail <- 0.5

lower_tail <- function(x, q = q_tail) {
  rowSums(x) <= q          # below red dashed line: u + v <= 0.5
}

upper_tail <- function(x, q = q_tail) {
  rowSums(x) >= 2 - q      # above blue dashed line: u + v >= 1.5
}

lt_id <- lower_tail(xtc_l[time, ])
ut_id <- upper_tail(xtc_r[time, ])

#===========================================================
# Optional: check correlations
#===========================================================

cor.test(f_cop[,1], f_cop[,2], method = "spearman")$estimate
cor.test(xtc_l[time,1], xtc_l[time,2], method = "spearman")$estimate
cor.test(xtc_r[time,1], xtc_r[time,2], method = "spearman")$estimate

#===========================================================
# Helper function: add tail threshold lines
#===========================================================

add_tail_lines <- function() {
  abline(a = 0.5, b = -1, col = "red",  lty = 2, lwd = 1)
  abline(a = 1.0, b = -1, col = "grey70", lty = 3, lwd = 1)
  abline(a = 1.5, b = -1, col = "blue", lty = 2, lwd = 1)
}

#===========================================================
# Full conceptual figure
#===========================================================

pdf(here("RESULTS/conceptual_fig_cop_ts.pdf"),
    height = 6, width = 6)

op <- par(mfrow = c(3, 2),
          mar = c(2.2, 2.2, 1.2, 1.2),
          mgp = c(1.4, 0.4, 0),
          pty = "s")

#-----------------------------------------------------------
# A1: no tail dependence
#-----------------------------------------------------------

plot(f_cop[,1], f_cop[,2],
     col = "black", bg = rgb(0, 0, 0, 0.3),
     pch = 21, cex = 1.4,
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n",
     ann = FALSE)

add_tail_lines()
text(0.1, 0.9, "A1", cex = 1.4)
text(0.58, 0.32, "No\ntail-dependence", cex = 0.8)

#-----------------------------------------------------------
# A2: no tail-dependence time series
#-----------------------------------------------------------

plot(1:N, f_cop[,1],
     col = "black", pch = 1,
     cex = 1.2, lwd = 0.6,
     xlim = c(1, N), ylim = c(0, 1),
     xaxt = "n", yaxt = "n",
     ann = FALSE, type = "b")

lines(1:N, f_cop[,2],
      col = rgb(0, 0, 0, 0.45),
      pch = 16, cex = 0.8, lwd = 0.6,
      type = "b")

text(35, 0.92, "A2", cex = 1.4)

#-----------------------------------------------------------
# B1: lower-tail dependence
# only lower-tail points are red
#-----------------------------------------------------------

plot(xtc_l[time,1], xtc_l[time,2],
     col = "grey50", bg = rgb(0, 0, 0, 0.12),
     pch = 21, cex = 1.4,
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n",
     ann = FALSE)

points(xtc_l[time,1][lt_id], xtc_l[time,2][lt_id],
       col = "red", bg = rgb(1, 0, 0, 0.3),
       pch = 21, cex = 1.4)

add_tail_lines()
text(0.82, 0.12, "B1", col = "red", cex = 1.4)
text(0.54, 0.18, "Lower\ntail-dependence", col = "red", cex = 0.8)

#-----------------------------------------------------------
# B2: lower-tail time series
# only lower-tail years are red
#-----------------------------------------------------------

plot(1:N, xtc_l[time,1],
     col = "grey70", pch = 1,
     cex = 1.2, lwd = 0.6,
     xlim = c(1, N), ylim = c(0, 1),
     xaxt = "n", yaxt = "n",
     ann = FALSE, type = "b")

lines(1:N, xtc_l[time,2],
      col = rgb(0, 0, 0, 0.25),
      pch = 16, cex = 0.8, lwd = 0.6,
      type = "b")

points((1:N)[lt_id], xtc_l[time,1][lt_id],
       col = "red", pch = 1, cex = 1.2)

points((1:N)[lt_id], xtc_l[time,2][lt_id],
       col = "red", pch = 16, cex = 0.8)

text(35, 0.92, "B2", col = "red", cex = 1.4)

#-----------------------------------------------------------
# C1: upper-tail dependence
# only upper-tail points are blue
#-----------------------------------------------------------

plot(xtc_r[time,1], xtc_r[time,2],
     col = "grey50", bg = rgb(0, 0, 0, 0.12),
     pch = 21, cex = 1.4,
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n",
     ann = FALSE)

points(xtc_r[time,1][ut_id], xtc_r[time,2][ut_id],
       col = "blue", bg = rgb(0, 0, 1, 0.3),
       pch = 21, cex = 1.4)

add_tail_lines()
text(0.82, 0.12, "C1", col = "blue", cex = 1.4)
text(0.28, 0.78, "Upper\ntail-dependence", col = "blue", cex = 0.8)

#-----------------------------------------------------------
# C2: upper-tail time series
# only upper-tail years are blue
#-----------------------------------------------------------

plot(1:N, xtc_r[time,1],
     col = "grey70", pch = 1,
     cex = 1.2, lwd = 0.6,
     xlim = c(1, N), ylim = c(0, 1),
     xaxt = "n", yaxt = "n",
     ann = FALSE, type = "b")

lines(1:N, xtc_r[time,2],
      col = rgb(0, 0, 0, 0.25),
      pch = 16, cex = 0.8, lwd = 0.6,
      type = "b")

points((1:N)[ut_id], xtc_r[time,1][ut_id],
       col = "blue", pch = 1, cex = 1.2)

points((1:N)[ut_id], xtc_r[time,2][ut_id],
       col = "blue", pch = 16, cex = 0.8)

text(35, 0.92, "C2", col = "blue", cex = 1.4)

par(op)
dev.off()