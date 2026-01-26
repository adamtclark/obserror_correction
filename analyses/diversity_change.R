# The following code replicates Fig. 5 in the manuscript

#setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
#source("util/save_model_outputs.R")
require(brms)
require(nlme)
load("output/brms_coverid_models_small.rda")
#load("output/diversity_change.rda")

set.seed(123432)

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.5
cex_level = 0.8
lwduse = 1.5

source("util/diversity_change_functions.R")

# loop: test observed change across deltaS values
deltaS_lvls = -10:10

simout_rich = NULL
simout_shannon = NULL
simout_simpson = NULL
for(i in 1:length(deltaS_lvls)) {
  divdat_new = simulate_change(divdat, deltaS = deltaS_lvls[i])
  
  divdat_0 = simulate_noise(divdat)
  divdat_1 = simulate_noise(divdat_new)
  
  beta_out_rich = get_beta(divdat_out_1 = divdat_0, divdat_out_2 = divdat_1, q = 0)
  simout_rich = rbind(simout_rich,
                 data.frame(beta_out_rich, dS = deltaS_lvls[i]))
  
  beta_out_shannon = get_beta(divdat_out_1 = divdat_0, divdat_out_2 = divdat_1, q = 1)
  simout_shannon = rbind(simout_shannon,
                      data.frame(beta_out_shannon, dS = deltaS_lvls[i]))
  
  beta_out_simpson = get_beta(divdat_out_1 = divdat_0, divdat_out_2 = divdat_1, q = 2)
  simout_simpson = rbind(simout_simpson,
                      data.frame(beta_out_simpson, dS = deltaS_lvls[i]))
  
  cat("\r", round(i/length(deltaS_lvls),2))
}

## analyse and plot results
plot_delta_function = function(simout, dobs, colu = 1, ylim=c(-1,1), doplot=TRUE) {
  simout$derror = dobs
  mod = lme(derror~-1+as.factor(dS),
            data = simout,
            random = ~1|SITE_CODE/plot)
  cfs = summary(mod)$tTable
  dS_lst = sort(unique(simout$dS))
  
  if(doplot) {
    plot(0, 0,
         type = "n",,
         xlim = range(dS_lst),
         ylim = ylim,
         ylab = "", xlab = "")
    abline(h=0, lty=2)
    abline(v=0,lty=3)
  }
  polygon(c(dS_lst, rev(dS_lst)),
          c(cfs[,1]+cfs[,2]*qnorm(0.025),
            rev(cfs[,1]+cfs[,2]*qnorm(0.975))),
          col = adjustcolor(colu,alpha.f = alpha_level),
          border = NA)
  lines(dS_lst, cfs[,1], lwd = lwduse,
        col = adjustcolor(colu,alpha.f = alpha_level))
}

annotate_plot = function(xlm = c(-10, 10), ylm = c(-0.14, 0.14)) {
  text(xlm[1], ylm[2], "underestimate species loss", cex = 0.8, pos = 4)
  text(xlm[2], ylm[1], "underestimate species gain", cex = 0.8, pos = 2)
  text(xlm[1], ylm[1], "overestimate species loss", cex = 0.8, pos = 4)
  text(xlm[2], ylm[2], "overestimate species gain", cex = 0.8, pos = 2)
}

annotate_plot2 = function(xlm = c(-10, 10), ylm = c(-0.04, 0.14)) {
  text(xlm[2], ylm[2], "overestimate turnover", cex = 0.8, pos = 2)
  text(xlm[2], ylm[1], "underestimate turnover", cex = 0.8, pos = 2)
}

pdf("figures/diversity_change.pdf", width = 7, height = 6)
#png("figures/diversity_change.png", width = 7, height = 6, units = "in", res = 200)

par(mfrow = c(3,2), mar=c(2,2,1.5,1), oma = c(2,2,0,0))
lm_alpha = c(-0.15, 0.15)
lm_beta = c(-0.05, 0.15)

#TOTAL:
#rich
dtrue = (simout_rich$alpha_2_true-simout_rich$alpha_1_true)
dobs_same = ((simout_rich$alpha_2_same-simout_rich$alpha_1_same)-dtrue)/simout_rich$alpha_1_same
plot_delta_function(simout = simout_rich,dobs = dobs_same, ylim = lm_alpha, colu = collst[3])
annotate_plot()

dobs_diff = ((simout_rich$alpha_2_diff-simout_rich$alpha_1_diff)-dtrue)/simout_rich$alpha_1_diff
plot_delta_function(simout = simout_rich,dobs = dobs_diff,ylim=lm_alpha, doplot = FALSE, colu = collst[4])
title("A. Richness, alpha", adj = 0)

dtrue = (simout_rich$beta_true)
dobs_same = (simout_rich$beta_same)-dtrue
plot_delta_function(simout = simout_rich,dobs = dobs_same, ylim = lm_beta, colu = collst[3])
annotate_plot2()

dobs_diff = (simout_rich$beta_different)-dtrue
plot_delta_function(simout = simout_rich,dobs = dobs_diff, ylim = lm_beta, colu = collst[4], doplot = FALSE)
title("D. Richness, beta", adj = 0)

# shannon
dtrue = (simout_shannon$alpha_2_true-simout_shannon$alpha_1_true)
dobs_same = ((simout_shannon$alpha_2_same-simout_shannon$alpha_1_same)-dtrue)/simout_shannon$alpha_1_same
plot_delta_function(simout = simout_shannon,dobs = dobs_same, ylim = lm_alpha, colu = collst[3])
annotate_plot()

dobs_diff = ((simout_shannon$alpha_2_diff-simout_shannon$alpha_1_diff)-dtrue)/simout_shannon$alpha_1_diff
plot_delta_function(simout = simout_shannon,dobs = dobs_diff,ylim=lm_alpha, colu = collst[4], doplot = FALSE)
title("B. Shannon, alpha", adj = 0)

dtrue = (simout_shannon$beta_true)
dobs_same = (simout_shannon$beta_same)-dtrue
plot_delta_function(simout = simout_shannon,dobs = dobs_same, ylim = lm_beta, colu = collst[3])
annotate_plot2()

dobs_diff = (simout_shannon$beta_different)-dtrue
plot_delta_function(simout = simout_shannon,dobs = dobs_diff, ylim = lm_beta, colu = collst[4], doplot = FALSE)
title("E. Shannon, beta", adj = 0)

# simpson
dtrue = (simout_simpson$alpha_2_true-simout_simpson$alpha_1_true)
dobs_same = ((simout_simpson$alpha_2_same-simout_simpson$alpha_1_same)-dtrue)/simout_simpson$alpha_1_same
plot_delta_function(simout = simout_simpson,dobs = dobs_same, ylim = lm_alpha, colu = collst[3])
annotate_plot()

dobs_diff = ((simout_simpson$alpha_2_diff-simout_simpson$alpha_1_diff)-dtrue)/simout_simpson$alpha_1_diff
plot_delta_function(simout = simout_simpson,dobs = dobs_diff,ylim=lm_alpha, colu = collst[4], doplot = FALSE)
title("C. Simpson, alpha", adj = 0)

dtrue = (simout_simpson$beta_true)
dobs_same = (simout_simpson$beta_same)-dtrue
plot_delta_function(simout = simout_simpson,dobs = dobs_same, ylim = lm_beta, colu = collst[3])
annotate_plot2()

dobs_diff = (simout_simpson$beta_different)-dtrue
plot_delta_function(simout = simout_simpson,dobs = dobs_diff, ylim = lm_beta, colu = collst[4], doplot = FALSE)
title("F. Simpson, beta", adj = 0)

mtext("Relative bias", side = 2, outer = TRUE, line = 0.6)
mtext("Simulated change in richness", side = 1, outer = TRUE, line = 0.6)

legend(-10, .16, c("Same Surveyor", "Different Surveyors"),
       #fill = adjustcolor(collst[3:4], alpha_level),
       lty = 1, col = collst[3:4],
       bty = "n", cex = 1.2)

dev.off()

