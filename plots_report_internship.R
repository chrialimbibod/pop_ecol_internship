###############################################################################
# Script to visualize results from wild boar simulation
#
# Christoph Imboden
#
# Created: 26.06.2024 CAI
#
###############################################################################

# clear R's memory :D
rm(list=ls())

# Load libraries that might be needed
library(viridis) # must in this script
library(scales)
library(RColorBrewer)
library(knitr) # only needed when doing stuff with matrices etc. & LaTex
library(plyr) # must in this script

# set the working directory (if required/wanted)
# setwd("~/.../.../")

# if utils are downloaded -> source("~/.../.../")
# convenient to "load" the function

##### 1) Set parameters #####

# t_span used
t_span <- 25

# N0 used
N0 <- rep(100, 6)

## -> became in the simulations:
## [1] 107 113 105 110  95  98

##### 2) Hunting regime 1 #####
# Simulate
set.seed(00072024) # set seed for reproducibility (does not apply to random realization!)
reg1 <- sim_population_hunted_JAGS(N0 = N0, t_span = t_span, hunting_regime = 1)

# Plot for random realization annual lambda and mean lambda
ylim <- c(0.5, 1.5)
plot(x = 1:t_span, y = reg1$JAGS_output$sims.list$lambda.annual[reg1$rand.realization,], type = "l", xlab = "Years", ylab = "",
     axes = F, ylim = ylim)
abline(h = reg1$JAGS_output$mean$lambda, lwd = 2, col = "red")
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)
legend(x= 1, y = 1.5, legend = c("Annual lambda", "Mean lambda"), lty = 1, col = c("black", "red"), bty = "n", lwd = c(1, 1), cex = 0.8)

# Plot Ntot of first 500 realizations 
colos <- turbo(499)
ylim <- c(0, round_any(max(reg1$JAGS_output$sims.list$Ntot[1:500,]), 1000, ceiling))
plot(x = 1:(t_span+1), y = reg1$JAGS_output$sims.list$Ntot[1,], type = "l", ylim = ylim,
     xlab = "Years", ylab = "", axes = F, main = "Total abundance in first 500 realizations", cex.main = 1)
for (i in 2:500) {
  lines(x = 1:(t_span+1), y = reg1$JAGS_output$sims.list$Ntot[i,], type = "l", col = colos[i])
}
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'") , las=1)

# Plot for N available prior to hunting, N hunted and N alive
# "avail" can be calculated with function "calc.avail" (available in utils)
ylim <- c(0, round_any(max(reg1$N.tot.t.random.realization, reg1$N.hunted.t.random.realization, avail, na.rm = T), 1000, ceiling))
plot(x = 1:(t_span), y = reg1$N.tot.t.random.realization[-1], type = "b", pch = 16, ylim = ylim,
     xlab = "Years", ylab = "", axes = F)
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, at = pretty(ylim), labels = format(pretty(0:max(avail)), big.mark = "'"), las=1)
lines(x = (1:t_span)-0.5, y = avail, type = "b", pch = 16, col = "blue")
lines(x = (1:t_span)-0.25, y = c(reg1$N.hunted.t.random.realization), type = "b", pch = 16, col = "red")
legend(x = 14.5, y = ylim[2], c("Total number of individuals","Total number available prior to hunting", 'Total number hunted'), pch=c(16, 16, 16),
       col=c('black','blue', 'red'), bty='n', lwd=c(1, 1, 1), cex = 0.8)

# Plot for N hunted animals per sex- and weight-class
ylim <- c(0, max(reg1$HR_df.random.realization))
colos <- turbo(6)
plot(x = 1:t_span, reg1$HR_df.random.realization[1,], type = "b", pch = 16, ylim = ylim, xlab = "Years", ylab = "", axes = F,
     main = paste("Number of hunted per sex and size class in realization", reg1$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)
points(x = 1:t_span, reg1$HR_df.random.realization[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:t_span, reg1$HR_df.random.realization[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:t_span, reg1$HR_df.random.realization[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:t_span, reg1$HR_df.random.realization[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:t_span, reg1$HR_df.random.realization[6,], type = "b", pch = 16, col = colos[6])
legend('bottom', c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

# Plot for mean N hunted animals per sex- and weight-class
ylim <- c(0, 750)
colos <- turbo(6)
plot(x = 1:t_span, reg1$mean_HR_df[1,], type = "b", pch = 16, ylim = ylim, xlab = "Years", ylab = "", axes = F,
     main = "Mean number of hunted per sex and size class", cex.main = 1, col = colos[1])
segments(1:t_span, reg1$JAGS_output$q2.5$hunt_f_s, 1:t_span, reg1$JAGS_output$q97.5$hunt_f_s, col=colos[1])
points(x = (1:t_span)+0.1, reg1$mean_HR_df[2,], type = "b", pch = 16, col = colos[2])
segments((1:t_span)+0.1, reg1$JAGS_output$q2.5$hunt_f_m, (1:t_span)+0.1, reg1$JAGS_output$q97.5$hunt_f_m, col=colos[2])
points(x = (1:t_span)-0.1, reg1$mean_HR_df[3,], type = "b", pch = 16, col = colos[3])
segments((1:t_span)-0.1, reg1$JAGS_output$q2.5$hunt_f_l, (1:t_span)-0.1, reg1$JAGS_output$q97.5$hunt_f_l, col=colos[3])
points(x = (1:t_span)-0.25, reg1$mean_HR_df[4,], type = "b", pch = 16, col = colos[4])
segments((1:t_span)-0.25, reg1$JAGS_output$q2.5$hunt_m_s, (1:t_span)-0.25, reg1$JAGS_output$q97.5$hunt_m_s, col=colos[4])
points(x = (1:t_span)+0.25, reg1$mean_HR_df[5,], type = "b", pch = 16, col = colos[5])
segments((1:t_span)+0.25, reg1$JAGS_output$q2.5$hunt_m_m, (1:t_span)+0.25, reg1$JAGS_output$q97.5$hunt_m_m, col=colos[5])
points(x = (1:t_span)+0.4, reg1$mean_HR_df[6,], type = "b", pch = 16, col = colos[6])
segments((1:t_span)+0.4, reg1$JAGS_output$q2.5$hunt_m_l, (1:t_span)+0.4, reg1$JAGS_output$q97.5$hunt_m_l, col=colos[6])
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)

##### 3) Hunting regime 2 #####
# Simulate
set.seed(00072024)
system.time(reg2 <- sim_population_hunted_JAGS(N0 = N0, t_span = t_span, hunting_regime = 2))

# Plot Ntot of first 500 realizations 
colos <- turbo(499)
ylim <- c(0, round_any(max(reg2$JAGS_output$sims.list$Ntot[1:500,]), 10000, ceiling))
plot(x = 1:(t_span+1), y = reg2$JAGS_output$sims.list$Ntot[1,], type = "l", ylim = ylim,
     xlab = "Years", ylab = "", axes = F, main = "Total abundance in first 500 realizations", cex.main = 1)
for (i in 2:500) {
  lines(x = 1:(t_span+1), y = reg2$JAGS_output$sims.list$Ntot[i,], type = "l", col = colos[i])
}
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'") , las=1)

# Plot for N available prior to hunting, N hunted and N alive
# "avail" can be calculated with function "calc.avail" (available in utils)
ylim <- c(0, round_any(max(reg2$N.tot.t.random.realization, reg2$N.hunted.t.random.realization, avail, na.rm = T), 1000, ceiling))
plot(x = 1:(t_span), y = reg2$N.tot.t.random.realization[-1], type = "b", pch = 16, ylim = ylim,
     xlab = "Years", ylab = "", axes = F)
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'"), las=1)
lines(x = (1:(t_span))-0.5, y = avail, type = "b", pch = 16, col = "blue")
lines(x = (1:(t_span))-0.25, y = c(reg2$N.hunted.t.random.realization), type = "b", pch = 16, col = "red")
legend(x = 2, y = ylim[2], c("Total number of individuals", 'Total number available prior to hunting', 'Total number hunted'), pch=c(16, 16, 16),
       col=c('black', 'blue' ,'red'), bty='n', lwd=c(1, 1, 1), cex = 0.8)

# Plot for N hunted animals per sex- and weight-class
ylim <- c(0, round_any(max(reg2$HR_df.random.realization), 1000))
colos <- turbo(6)
plot(x = 1:t_span, reg2$HR_df.random.realization[1,], type = "b", pch = 16, ylim = ylim, xlab = "Years", ylab = "", axes = F,
     main = paste("Number of hunted per sex and size class in realization", reg2$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(from = 0, to = ylim[2], by = 100), tcl = -0.25, labels = NA)
axis(2, at = pretty(0:ylim[2]), labels = format(pretty(0:ylim[2]), big.mark = "'"), las=1)
points(x = 1:t_span, reg2$HR_df.random.realization[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:t_span, reg2$HR_df.random.realization[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:t_span, reg2$HR_df.random.realization[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:t_span, reg2$HR_df.random.realization[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:t_span, reg2$HR_df.random.realization[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

# Plot for N alive animals per sex- and weight-class
# Produce a dataframe
dd <- matrix(NA, nrow = 6, ncol = (t_span+1))
for (i in 1:(t_span+1)){
  dd[,i] <- reg2$JAGS_output$sims.list$N[reg2$rand.realization, ,i]
}

ylim <- c(0, round_any(max(dd), 1000))
colos <- turbo(6)
plot(x = 1:(t_span+1), dd[1,], type = "b", pch = 16, ylim = ylim, xlab = "Years", ylab = "", axes = F,
     main = paste("Number of animals per sex and size class in realization", reg2$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at = seq(from = 0, to = ylim[2], 250), tcl=-0.25, labels = NA)
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'") , las=1)
points(x = 1:(t_span+1), dd[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:(t_span+1), dd[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:(t_span+1), dd[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:(t_span+1), dd[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:(t_span+1), dd[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = (max(dd)+200), c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

# Visualize the stage distribution (not in report)
r <- reg2$rand.realization
# Produce a dataframe
dssd <- matrix(NA, nrow = 6, ncol = t_span)
for (t in 1:t_span) {
  dssd[,t] <- reg2$JAGS_output$sims.list$stage.distr[r, ,t]
}
ylim <- c(0,1)
plot(x = 1:25, y = dssd[1,], type = "b", pch = 16, col = colos[1], axes = F, xlab = "Years", ylab = "", cex = 0.8, ylim = ylim,
     main = "Proportion sex- and weight-class of total population", cex.main = 1)
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, at = seq(from = 0, to = ylim[2], 0.1), tcl=-0.25, labels = NA)
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'") , las=1)
points(1:25, y = dssd[2,], type = "b", pch = 16, col = colos[2], cex = 0.8)
points(1:25, y = dssd[3,], type = "b", pch = 16, col = colos[3], cex = 0.8)
points(1:25, y = dssd[4,], type = "b", pch = 16, col = colos[4], cex = 0.8)
points(1:25, y = dssd[5,], type = "b", pch = 16, col = colos[5], cex = 0.8)
points(1:25, y = dssd[6,], type = "b", pch = 16, col = colos[6], cex = 0.8)
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

##### 4) Hunting regime 3 #####
# Simulate
set.seed(00072024)
system.time(reg3 <- sim_population_hunted_JAGS(N0 = N0, t_span = t_span, hunting_regime = 3))

# Plot Ntot of first 500 realizations
colos <- turbo(499)
ylim <- c(0, round_any(max(reg3$JAGS_output$sims.list$Ntot[1:500,]), 1000))
plot(x = 1:(t_span+1), y = reg3$JAGS_output$sims.list$Ntot[1,], type = "l", ylim = ylim,
     xlab = "Years", ylab = "", axes = F, main = "Total abundance in first 500 realizations", cex.main = 1)
for (i in 2:500) {
  lines(x = 1:(t_span+1), y = reg3$JAGS_output$sims.list$Ntot[i,], type = "l", col = colos[i])
}
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, at = seq(from = 0, to = ylim[2], 500), tcl=-0.25, labels = NA)
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'") , las=1)

# Plot for N available prior to hunting, N hunted and N alive
# "avail" can be calculated with function "calc.avail" (available in utils)
ylim <- c(0, round_any(max(reg3$N.tot.t.random.realization, reg3$N.hunted.t.random.realization, avail, na.rm = T), 10000))
plot(x = 1:(t_span), y = reg3$N.tot.t.random.realization[-1], type = "b", pch = 16, ylim = ylim,
     xlab = "Years", ylab = "", axes = F)
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at = pretty(avail), labels = pretty(avail), las=1)
axis(2, at = pretty(avail), labels = format(pretty(avail), big.mark = "'"), las=1)
lines(x = (1:(t_span))-0.5, y = avail, type = "b", pch = 16, col = "blue")
lines(x = (1:(t_span))-0.25, y = c(reg3$N.hunted.t.random.realization), type = "b", pch = 16, col = "red")
legend(x = 2, y = ylim[2], c("Total number of individuals", 'Total number available prior to hunting', 'Total number hunted'), pch=c(16, 16, 16),
       col=c('black', 'blue' ,'red'), bty='n', lwd=c(1, 1, 1), cex = 0.8)

# Plot for N hunted animals per sex- and weight-class
ylim <- c(0, round_any(max(reg3$HR_df.random.realization), 1000))
colos <- turbo(6)
plot(x = 1:t_span, reg3$HR_df.random.realization[1,], type = "b", pch = 16, ylim = (ylim), xlab = "Years", ylab = "", axes = F,
     main = paste("Number of hunted per sex and size class in realization", reg3$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
#axis(2, at=seq(0, ylim[2], 50), tcl=-0.25, labels=NA)
axis(2, at = pretty(ylim), tcl=-0.5, labels = format(pretty(ylim), big.mark = "'"), las = 1)
points(x = 1:t_span, reg3$HR_df.random.realization[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:t_span, reg3$HR_df.random.realization[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:t_span, reg3$HR_df.random.realization[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:t_span, reg3$HR_df.random.realization[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:t_span, reg3$HR_df.random.realization[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = (max(reg3$HR_df.random.realization)+100), c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

# Plot for N alive animals per sex- and weight-class
# Produce a dataframe
ddr3 <- matrix(NA, nrow = 6, ncol = (t_span+1))
for (i in 1:(t_span+1)){
  ddr3[,i] <- reg3$JAGS_output$sims.list$N[reg3$rand.realization, ,i]
}


ylim <- c(0, round_any(max(ddr3), 1000, ceiling))
colos <- turbo(6)
plot(x = 1:(t_span+1), ddr3[1,], type = "b", pch = 16, ylim = ylim, xlab = "Years", ylab = "", axes = F,
     main = paste("Number of animals per sex and size class in realization", reg3$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at = seq(0, ylim[2], 100), tcl=-0.25, labels=NA)
axis(2, at = pretty(ylim), tcl=-0.5, labels=format(pretty(ylim), big.mark = "'"), las=1)
points(x = 1:(t_span+1), ddr3[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:(t_span+1), ddr3[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:(t_span+1), ddr3[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:(t_span+1), ddr3[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:(t_span+1), ddr3[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

##### 5) Hunting regime 4 #####
# Simulate
set.seed(00072024)
system.time(reg4 <- sim_population_hunted_JAGS(N0 = N0, t_span = t_span, hunting_regime = 4))

# Plot Ntot of first 500 realizations
colos <- turbo(499)
ylim <- c(0, round_any(max(reg4$JAGS_output$sims.list$Ntot[1:500,]), 1000))
plot(x = 1:(t_span+1), y = reg4$JAGS_output$sims.list$Ntot[1,], type = "l", ylim = ylim,
     xlab = "Years", ylab = "", axes = F, main = "Total abundance in first 500 realizations", cex.main = 1)
for (i in 2:500) {
  lines(x = 1:(t_span+1), y = reg4$JAGS_output$sims.list$Ntot[i,], type = "l", col = colos[i])
}
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(0, ylim[2], 250), tcl=-0.25, labels=NA)
axis(2, at=pretty(ylim), tcl=-0.5, labels=format(pretty(ylim), big.mark = "'"), las=1)
# abline(h=sum(reg4$JAGS_output$sims.list$N[1, ,1]), lwd = 3)
# lowers <- which(reg4$JAGS_output$sims.list$lambda < 1)
# lowers_colos <- viridis(length(lowers))
# c <- 1
# for (i in lowers) {
#   lines(x = 1:(t_span+1), y = reg4$JAGS_output$sims.list$Ntot[i,], type = "l", col = lowers_colos[c])
#   c <- c + 1
# }

# Plot for N available prior to hunting, N hunted and N alive
# "avail" can be calculated with function "calc.avail" (available in utils)
ylim <- c(0, round_any(max(reg4$N.tot.t.random.realization, reg4$N.hunted.t.random.realization, avail, na.rm = T), 1000, ceiling))
plot(x = 1:(t_span), y = reg4$N.tot.t.random.realization[-1], type = "b", pch = 16, ylim = ylim,
     xlab = "Years", ylab = "", axes = F)
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at = pretty(avail), labels = pretty(avail), las=1)
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'"), las=1)
lines(x = (1:(t_span))-0.5, y = avail, type = "b", pch = 16, col = "blue")
lines(x = (1:(t_span))-0.25, y = c(reg4$N.hunted.t.random.realization), type = "b", pch = 16, col = "red")
legend(x = 2, y = ylim[2], c("Total number of individuals", 'Total number available prior to hunting', 'Total number hunted'), pch=c(16, 16, 16),
       col=c('black', 'blue' ,'red'), bty='n', lwd=c(1, 1, 1), cex = 0.8)

# Plot for N hunted animals per sex- and weight-class
ylim <- c(0, round_any(max(reg4$HR_df.random.realization), 100))
colos <- turbo(6)
plot(x = 1:t_span, reg4$HR_df.random.realization[1,], type = "b", pch = 16, ylim = (ylim), xlab = "Years", ylab = "", axes = F,
     main = paste("Number of hunted per sex and size class in realization", reg4$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(0, ylim[2], 20), tcl=-0.25, labels=NA)
axis(2, at = pretty(ylim), tcl=-0.5, labels=format(pretty(ylim), big.mark = "'"), las=1)
points(x = 1:t_span, reg4$HR_df.random.realization[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:t_span, reg4$HR_df.random.realization[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:t_span, reg4$HR_df.random.realization[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:t_span, reg4$HR_df.random.realization[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:t_span, reg4$HR_df.random.realization[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

# Plot for N alive animals per sex- and weight-class
# Produce a dataframe
ddr4 <- matrix(NA, nrow = 6, ncol = (t_span+1))
for (i in 1:(t_span+1)){
  ddr4[,i] <- reg4$JAGS_output$sims.list$N[reg4$rand.realization, ,i]
}

ylim <- c(0, round_any(max(ddr4), 1000))
colos <- turbo(6)
plot(x = 1:(t_span+1), ddr4[1,], type = "b", pch = 16, ylim = ylim, xlab = "Years", ylab = "", axes = F,
     main = paste("Number of animals per sex and size class in realization", reg4$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, at = pretty(ylim), tcl=-0.5, labels=format(pretty(ylim), big.mark = "'"), las=1)
points(x = 1:(t_span+1), ddr4[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:(t_span+1), ddr4[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:(t_span+1), ddr4[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:(t_span+1), ddr4[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:(t_span+1), ddr4[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

##### 6) Hunting regime 5 #####
# Simulate
set.seed(00072024)
system.time(reg5 <- sim_population_hunted_JAGS(N0 = N0, t_span = t_span, hunting_regime = 5))

# Plot Ntot of first 500 realizations
colos <- turbo(499)
ylim <- c(0, round_any(max(reg5$JAGS_output$sims.list$Ntot[1:500,]), 1000, floor))
plot(x = 1:(t_span+1), y = reg5$JAGS_output$sims.list$Ntot[1,], type = "l", ylim = ylim,
     xlab = "Years", ylab = "", axes = F, main = "Total abundance in first 500 realizations", cex.main = 1)
for (i in 2:500) {
  lines(x = 1:(t_span+1), y = reg5$JAGS_output$sims.list$Ntot[i,], type = "l", col = colos[i])
}
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(0, ylim[2], 500), tcl=-0.25, labels=NA)
axis(2, at=pretty(ylim), tcl=-0.5, labels=format(pretty(ylim), big.mark = "'"), las=1)

# Plot for N available prior to hunting, N hunted and N alive
# "avail" can be calculated with function "calc.avail" (available in utils)
ylim <- c(0, round_any(max(reg5$N.tot.t.random.realization, reg5$N.hunted.t.random.realization, avail, na.rm = T), 1000, ceiling))
plot(x = 1:(t_span), y = reg5$N.tot.t.random.realization[-1], type = "b", pch = 16, ylim = ylim,
     xlab = "Years", ylab = "", axes = F)
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at = pretty(avail), labels = pretty(avail), las=1)
# axis(2, at=seq(0, ylim[2], 500), tcl=-0.25, labels=NA)
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'"), tcl=-0.5, las=1)
lines(x = (1:(t_span))-0.5, y = avail, type = "b", pch = 16, col = "blue")
lines(x = (1:(t_span))-0.25, y = c(reg5$N.hunted.t.random.realization), type = "b", pch = 16, col = "red")
legend(x = 2, y = ylim[2], c("Total number of individuals", 'Total number available prior to hunting', 'Total number hunted'), pch=c(16, 16, 16),
       col=c('black', 'blue' ,'red'), bty='n', lwd=c(1, 1, 1), cex = 0.8)

# Plot for N hunted animals per sex- and weight-class
ylim <- c(0, round_any(max(reg5$HR_df.random.realization), 1500, ceiling))
colos <- turbo(6)
plot(x = 1:t_span, reg5$HR_df.random.realization[1,], type = "b", pch = 16, ylim = (ylim), xlab = "Years", ylab = "", axes = F,
     main = paste("Number of hunted per sex and size class in realization", reg5$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(0, ylim[2], 100), tcl=-0.25, labels=NA)
axis(2, at=pretty(ylim), labels = format(pretty(ylim), big.mark = "'"), tcl=-0.5, las=1)
points(x = 1:t_span, reg5$HR_df.random.realization[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:t_span, reg5$HR_df.random.realization[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:t_span, reg5$HR_df.random.realization[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:t_span, reg5$HR_df.random.realization[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:t_span, reg5$HR_df.random.realization[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)
# Plot for N alive animals per sex- and weight-class
# Produce a dataframe
ddr5 <- matrix(NA, nrow = 6, ncol = (t_span+1))
for (i in 1:(t_span+1)){
  ddr5[,i] <- reg5$JAGS_output$sims.list$N[reg5$rand.realization, ,i]
}

ylim <- c(0, round_any(max(ddr5), 1000))
colos <- turbo(6)
plot(x = 1:(t_span+1), ddr5[1,], type = "b", pch = 16, ylim = ylim, xlab = "Years", ylab = "", axes = F,
     main = paste("Number of animals per sex and size class in realization", reg5$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(0, ylim[2], 100), tcl=-0.25, labels=NA)
axis(2, at= pretty(ylim), tcl=-0.5, labels = format(pretty(ylim), big.mark = "'"), las=1)
points(x = 1:(t_span+1), ddr5[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:(t_span+1), ddr5[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:(t_span+1), ddr5[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:(t_span+1), ddr5[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:(t_span+1), ddr5[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

##### 7) Hunting regime 6 #####
# Simulate
set.seed(00072024)
system.time(reg6 <- sim_population_hunted_JAGS(N0 = N0, t_span = t_span, hunting_regime = 6))

# Plot Ntot of first 500 realizations
colos <- turbo(499)
ylim <- c(0, round_any(max(reg6$JAGS_output$sims.list$Ntot[1:500,]), 10000, ceiling))
plot(x = 1:(t_span+1), y = reg6$JAGS_output$sims.list$Ntot[1,], type = "l", ylim = ylim,
     xlab = "Years", ylab = "", axes = F, main = "Total abundance in first 500 realizations", cex.main = 1)
for (i in 2:500) {
  lines(x = 1:(t_span+1), y = reg6$JAGS_output$sims.list$Ntot[i,], type = "l", col = colos[i])
}
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(0, ylim[2], 1000), tcl=-0.25, labels=NA)
axis(2, at=pretty(ylim), labels = format(pretty(ylim), big.mark = "'"), tcl=-0.5, las=1)

# Plot for N available prior to hunting, N hunted and N alive
# "avail" can be calculated with function "calc.avail" (available in utils)
ylim <- c(0, round_any(max(reg6$N.tot.t.random.realization, reg6$N.hunted.t.random.realization, avail, na.rm = T), 1000, ceiling))
plot(x = 1:(t_span), y = reg6$N.tot.t.random.realization[-1], type = "b", pch = 16, ylim = ylim,
     xlab = "Years", ylab = "", axes = F)
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at = pretty(avail), labels = pretty(avail), las=1)
# axis(2, at=seq(0, ylim[2], 500), tcl=-0.25, labels=NA)
axis(2, at = pretty(ylim), labels = format(pretty(ylim), big.mark = "'"), tcl=-0.5, las=1)
lines(x = (1:(t_span))-0.5, y = avail, type = "b", pch = 16, col = "blue")
lines(x = (1:(t_span))-0.25, y = c(reg6$N.hunted.t.random.realization), type = "b", pch = 16, col = "red")
legend(x = 2, y = ylim[2], c("Total number of individuals", 'Total number available prior to hunting', 'Total number hunted'), pch=c(16, 16, 16),
       col=c('black', 'blue' ,'red'), bty='n', lwd=c(1, 1, 1), cex = 0.8)

# Plot for N hunted animals per sex- and weight-class
ylim <- c(0, round_any(max(reg6$HR_df.random.realization), 1000, ceiling))
colos <- turbo(6)
plot(x = 1:t_span, reg6$HR_df.random.realization[1,], type = "b", pch = 16, ylim = (ylim), xlab = "Years", ylab = "", axes = F,
     main = paste("Number of hunted per sex and size class in realization", reg6$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(0, ylim[2], 100), tcl=-0.25, labels=NA)
axis(2, at=pretty(ylim), tcl=-0.5, labels=format(pretty(ylim), big.mark = "'"), las=1)
points(x = 1:t_span, reg6$HR_df.random.realization[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:t_span, reg6$HR_df.random.realization[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:t_span, reg6$HR_df.random.realization[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:t_span, reg6$HR_df.random.realization[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:t_span, reg6$HR_df.random.realization[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

# Plot for N alive animals per sex- and weight-class
# Produce a dataframe
ddr6 <- matrix(NA, nrow = 6, ncol = (t_span+1))
for (i in 1:(t_span+1)){
  ddr6[,i] <- reg6$JAGS_output$sims.list$N[reg6$rand.realization, ,i]
}

ylim <- c(0, round_any(max(ddr6), 1000, ceiling))
colos <- turbo(6)
plot(x = 1:(t_span+1), ddr6[1,], type = "b", pch = 16, ylim = ylim, xlab = "Years", ylab = "", axes = F,
     main = paste("Number of animals per sex and size class in realization", reg6$rand.realization, sep = " "), cex.main = 1, col = colos[1])
axis(1, at=1:(t_span+1), tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
# axis(2, at=seq(0, ylim[2], 100), tcl=-0.25, labels=NA)
axis(2, at=pretty(ylim), labels = format(pretty(ylim), big.mark = "'"),  tcl=-0.5, las=1)
points(x = 1:(t_span+1), ddr6[2,], type = "b", pch = 16, col = colos[2])
points(x = 1:(t_span+1), ddr6[3,], type = "b", pch = 16, col = colos[3])
points(x = 1:(t_span+1), ddr6[4,], type = "b", pch = 16, col = colos[4])
points(x = 1:(t_span+1), ddr6[5,], type = "b", pch = 16, col = colos[5])
points(x = 1:(t_span+1), ddr6[6,], type = "b", pch = 16, col = colos[6])
legend(x = 1, y = ylim[2], c("Small females", "Medium females", "Large females", "Small males", "Medium males", "Large males"), pch=rep(16, 6),
       col=colos, bty='n', lwd=rep(1, 6), cex = 0.8, ncol = 2)

