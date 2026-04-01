# The following code replicates Fig. 1 in the manuscript

#setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
require(brms)
#load("output/cover_error.rda")

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.3
cex_level = 0.8
xsq = seq(0, 1, length = 1e3)

d = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
#d = d[d$SITE_CODE!="temple.us",]
d = d[d$trt=="Control",]

# get variance and mean
covdat = d[,grep("cover_", colnames(d))]
covdat[!is.na(covdat) & covdat==0] = NA
covdat = covdat/100
d$cov_var = apply(covdat, 1, function(x) var(x, na.rm=TRUE))
d$cov_mean = apply(covdat, 1, function(x) mean(x, na.rm=TRUE))
d$cov_cv = sqrt(d$cov_var)/d$cov_mean

mod_c0 <- brm(bf(cov_cv ~ cov_mean + (1|SITE_CODE),
                 hu ~ cov_mean + (1|SITE_CODE)),
              save_pars = save_pars(all = TRUE),
              data = d, family = hurdle_lognormal(), cores = 2)
mod_c0 # R-hats all < 1.01
#plot(mod_c0)
#pred_c = predict(mod_c0, re_formula = ~ 0, newdata = data.frame(cov_mean = xsq))
tmp = posterior_epred(mod_c0, re_formula = ~ 0, newdata = data.frame(cov_mean = xsq))
pred_c = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))))
colnames(pred_c) = c("Estimate", "Q2.5", "Q97.5")

par(mar=c(4,4,2,2))
plot((d$cov_mean), d$cov_cv,
     ylim = c(0, 2/sqrt(2)),
     xlim = c(0,1),
     cex = cex_level, col = adjustcolor(1, alpha_level),
     xlab = "Mean Cover", ylab = "Cover CV")
abline(h = 0, v = c(0,1), lty = 3)
#abline(a = 0, b = sqrt(2), lty = 3)
abline(h = 2/sqrt(2), lty=3)

matlines(xsq, pred_c[,c("Estimate", "Q2.5", "Q97.5")], lty = c(1,2,2), col = 1)

# graminoids vs. non-graminoids
mod_cg <- brm(bf(cov_cv ~ cov_mean*group + (1|SITE_CODE),
                 hu ~ cov_mean*group + (1|SITE_CODE)),
              save_pars = save_pars(all = TRUE),
              data = d, family = hurdle_lognormal(), cores = 2)
mod_cg # R-hats all <= 1.01
#plot(mod_cg)
#pred_cg = predict(mod_cg, re_formula = ~ 0, newdata = data.frame(cov_mean = xsq, group = "Graminoid"))
#pred_cng = predict(mod_cg, re_formula = ~ 0, newdata = data.frame(cov_mean = xsq, group = "Non-Graminoid"))

tmp = posterior_epred(mod_cg, re_formula = ~ 0, newdata = data.frame(cov_mean = xsq, group = "Graminoid"))
pred_cg = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))))
colnames(pred_cg) = c("Estimate", "Q2.5", "Q97.5")

tmp = posterior_epred(mod_cg, re_formula = ~ 0, newdata = data.frame(cov_mean = xsq, group = "Non-Graminoid"))
pred_cng = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))))
colnames(pred_cng) = c("Estimate", "Q2.5", "Q97.5")

par(mar=c(4,4,2,2))
plot((d$cov_mean), d$cov_cv,
     col = adjustcolor(collst[as.factor(d$group)], alpha.f = alpha_level),
     ylim = c(0, 2/sqrt(2)),
     xlim = c(0,1),
     cex = cex_level,
     xlab = "Mean Cover", ylab = "Cover CV")
abline(h = 0, v = c(0,1), lty = 3)
#abline(a = 0, b = sqrt(2), lty = 3)
abline(h = 2/sqrt(2), lty=3)

matlines(xsq, pred_cg[,c("Estimate", "Q2.5", "Q97.5")], lty = c(1,2,2), col = collst[1])
matlines(xsq, pred_cng[,c("Estimate", "Q2.5", "Q97.5")], lty = c(1,2,2), col = collst[2])

legend("topright", c("Graminoid", "Non-Graminoid"),
       lty = 1, pch = 1, col = collst[1:2], bty = "n")
loo(mod_c0, mod_cg)

# same surveyor
surveyor = d[,grep("surveyor_", colnames(d))]
ps = surveyor!=surveyor[,1]
ps[is.na(ps)] = TRUE

covdat_same = covdat
covdat_same[ps] = NA

d$cov_var_same = apply(covdat_same, 1, function(x) var(x, na.rm=TRUE))
d$cov_mean_same = apply(covdat_same, 1, function(x) mean(x, na.rm=TRUE))
d$cov_cv_same = sqrt(d$cov_var_same)/d$cov_mean_same

mod_c0s <- brm(bf(cov_cv_same ~ cov_mean_same + (1|SITE_CODE),
                  hu ~ cov_mean_same + (1|SITE_CODE)),
               save_pars = save_pars(all = TRUE),
               data = d, family = hurdle_lognormal(), cores = 2)
mod_c0s # R-hats all <= 1.01
#plot(mod_c0s)
#pred_cs = predict(mod_c0s, re_formula = ~ 0, newdata = data.frame(cov_mean_same = xsq), probs = c(0.025, pnorm(c(-1,1)), 0.975))
tmp = posterior_epred(mod_c0s, re_formula = ~ 0, newdata = data.frame(cov_mean_same = xsq))
pred_cs = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
colnames(pred_cs) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")


mod_cgs <- brm(bf(cov_cv_same ~ cov_mean_same*group + (1|SITE_CODE),
                  hu ~ cov_mean_same*group + (1|SITE_CODE)),
               save_pars = save_pars(all = TRUE),
               data = d, family = hurdle_lognormal(), cores = 2)
mod_cgs # R-hats all < 1.01
#plot(mod_cgs)
#pred_cgs = predict(mod_cgs, re_formula = ~ 0, newdata = data.frame(cov_mean_same = xsq, group = "Graminoid"), probs = c(0.025, pnorm(c(-1,1)), 0.975))
#pred_cngs = predict(mod_cgs, re_formula = ~ 0, newdata = data.frame(cov_mean_same = xsq, group = "Non-Graminoid"), probs = c(0.025, pnorm(c(-1,1)), 0.975))

tmp = posterior_epred(mod_cgs, re_formula = ~ 0, newdata = data.frame(cov_mean_same = xsq, group = "Graminoid"))
pred_cgs = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
colnames(pred_cgs) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_cgs, re_formula = ~ 0, newdata = data.frame(cov_mean_same = xsq, group = "Non-Graminoid"))
pred_cngs = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
colnames(pred_cngs) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")


# different surveyor
ps = surveyor==surveyor[,1]
ps[,-1][is.na(ps[,-1])] = TRUE
ps[,1][!is.na(ps[,1])] = FALSE

covdat_diff = covdat
covdat_diff[ps] = NA

d$cov_var_diff = apply(covdat_diff, 1, function(x) var(x, na.rm=TRUE))
d$cov_mean_diff = apply(covdat_diff, 1, function(x) mean(x, na.rm=TRUE))
d$cov_cv_diff = sqrt(d$cov_var_diff)/d$cov_mean_same

mod_c0d <- brm(bf(cov_cv_diff ~ cov_mean_diff + (1|SITE_CODE),
                  hu ~ cov_mean_diff + (1|SITE_CODE)),
               save_pars = save_pars(all = TRUE),
               data = d, family = hurdle_lognormal(), cores = 2)
mod_c0d # R-hats all <= 1.01
#plot(mod_c0d)
#pred_cd = predict(mod_c0d, re_formula = ~ 0, newdata = data.frame(cov_mean_diff = xsq), probs = c(0.025, pnorm(c(-1,1)), 0.975))
tmp = posterior_epred(mod_c0d, re_formula = ~ 0, newdata = data.frame(cov_mean_diff = xsq))
pred_cd = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
colnames(pred_cd) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

mod_cgd <- brm(bf(cov_cv_diff ~ cov_mean_diff*group + (1|SITE_CODE),
                  hu ~ cov_mean_diff*group + (1|SITE_CODE)),
               save_pars = save_pars(all = TRUE), 
               data = d, family = hurdle_lognormal(), cores = 2)
mod_cgd # R-hats all < 1.01
#plot(mod_cgd)
#pred_cgd = predict(mod_cgd, re_formula = ~ 0, newdata = data.frame(cov_mean_diff = xsq, group = "Graminoid"), probs = c(0.025, pnorm(c(-1,1)), 0.975))
#pred_cngd = predict(mod_cgd, re_formula = ~ 0, newdata = data.frame(cov_mean_diff = xsq, group = "Non-Graminoid"), probs = c(0.025, pnorm(c(-1,1)), 0.975))

tmp = posterior_epred(mod_cgd, re_formula = ~ 0, newdata = data.frame(cov_mean_diff = xsq, group = "Graminoid"))
pred_cgd = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
colnames(pred_cgd) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_cgd, re_formula = ~ 0, newdata = data.frame(cov_mean_diff = xsq, group = "Non-Graminoid"))
pred_cngd = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
colnames(pred_cngd) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")


# model comparison
loo_c0s = loo(mod_c0s, moment_match = TRUE) # Pareto k estimates < 0.7
loo_c0d = loo(mod_c0d, moment_match = TRUE) # Pareto k estimates < 0.7

loo_c0s_cgs = loo(mod_c0s, mod_cgs, moment_match = TRUE) # Pareto k estimates < 0.7; cgs best (-5.3, sd 5.6)
loo_c0d_cgd = loo(mod_c0d, mod_cgd, moment_match = TRUE) # Pareto k estimates < 0.7; cgd best (-13.5, sd 9.5)

# save models
#save(list = ls(), file = "output/cover_error.rda")

#####################
# plot total
pdf("figures/cover_error.pdf", width = 8, height = 4)
#png("figures/cover_error.png", width = 8, height = 4, units = "in", res = 200)

X = cbind(c(1,1), c(1,1), c(1,1), c(2,3),  c(2,3))
layout(X)
par(mar=c(2,2,2,2), oma = c(2,2,0,0))
plot((d$cov_mean_same), d$cov_cv_same,
     xaxs = "i",
     ylim = c(0, 2/sqrt(2)),
     xlim = c(0,1),
     cex = cex_level, col = adjustcolor(collst[3], alpha_level),
     xlab = "", ylab = "")
mtext("Mean non-zero species-level cover", 1, line = 2.4)
mtext("CV of non-zero species-level cover", 2, line = 2.4)
abline(h = 0, v = c(0,1), lty = 3)
#abline(a = 0, b = sqrt(2), lty = 3)
abline(h = 2/sqrt(2), lty=3)
lines(xsq, sqrt(0.01^2/2)/xsq, lty = 2, lwd = 1.5)

polygon(c(xsq, rev(xsq)), c(pred_cs[,c(3)], rev(pred_cs[,c(4)])), col = adjustcolor(collst[3], 0.5), border = NA)
lines(xsq, pred_cs[,1], lty = c(1), col = collst[3], lwd = 1.8)
#matlines(xsq, pred_cs[,c(1,4,5)], lty = c(1,2,2), col = collst[3])

points((d$cov_mean_diff), d$cov_cv_diff,
       pch = 2,
       cex = cex_level, col = adjustcolor(collst[4], alpha_level))
polygon(c(xsq, rev(xsq)), c(pred_cd[,c(3)], rev(pred_cd[,c(4)])), col = adjustcolor(collst[4], 0.5), border = NA)
lines(xsq, pred_cd[,1], lty = c(1), col = collst[4], lwd = 1.8)
#matlines(xsq, pred_cd[,c(1,4,5)], lty = c(1,2,2), col = collst[4])

legend(0.6, 1.4, c("Same Surveyor", "Different Surveyors"),
       lty = 1, pch = c(1,2), col = collst[3:4], bty = "n", cex = 1.2)
title("A. All Species", adj = 0)

# plot graminoids
ps = d$group == "Graminoid"
plot((d$cov_mean_same)[ps], d$cov_cv_same[ps],
     xaxs = "i",
     ylim = c(0, 2/sqrt(2)),
     xlim = c(0,1),
     cex = cex_level, col = adjustcolor(collst[3], alpha_level),
     xlab = "", ylab = "")
#mtext("Mean species-level cover", 1, line = 2.4)
#mtext("CV of species-level cover", 2, line = 2.4)
abline(h = 0, v = c(0,1), lty = 3)
#abline(a = 0, b = sqrt(2), lty = 3)
abline(h = 2/sqrt(2), lty=3)
lines(xsq, sqrt(0.01^2/2)/xsq, lty = 2, lwd = 1.5)
polygon(c(xsq, rev(xsq)), c(pred_cgs[,c(3)], rev(pred_cgs[,c(4)])), col = adjustcolor(collst[3], 0.5), border = NA)
lines(xsq, pred_cgs[,1], lty = c(1), col = collst[3], lwd = 1.8)
#matlines(xsq, pred_cgs[,c(1,4,5)], lty = c(1,2,2), col = collst[3])

points((d$cov_mean_diff)[ps], d$cov_cv_diff[ps],
       pch = 2,
       cex = cex_level, col = adjustcolor(collst[4], alpha_level))
polygon(c(xsq, rev(xsq)), c(pred_cgd[,c(3)], rev(pred_cgd[,c(4)])), col = adjustcolor(collst[4], 0.5), border = NA)
lines(xsq, pred_cgd[,1], lty = c(1), col = collst[4], lwd = 1.8)
#matlines(xsq, pred_cgd[,c(1,4,5)], lty = c(1,2,2), col = collst[4])

#legend("topright", c("Same Surveyor", "Different Surveyors"),
#       lty = 1, pch = c(1,2), col = collst[3:4], bty = "n")
title("B. Graminoids", adj = 0)


# plot non-graminoids
ps = d$group == "Non-Graminoid"
plot((d$cov_mean_same)[ps], d$cov_cv_same[ps],
     xaxs = "i",
     ylim = c(0, 2/sqrt(2)),
     xlim = c(0,1),
     cex = cex_level, col = adjustcolor(collst[3], alpha_level),
     xlab = "", ylab = "")

abline(h = 0, v = c(0,1), lty = 3)
#abline(a = 0, b = sqrt(2), lty = 3)
abline(h = 2/sqrt(2), lty=3)
lines(xsq, sqrt(0.01^2/2)/xsq, lty = 2, lwd = 1.5)

polygon(c(xsq, rev(xsq)), c(pred_cngs[,c(3)], rev(pred_cngs[,c(4)])), col = adjustcolor(collst[3], 0.5), border = NA)
lines(xsq, pred_cngs[,1], lty = c(1), col = collst[3], lwd = 1.8)
#matlines(xsq, pred_cngs[,c(1,4,5)], lty = c(1,2,2), col = collst[3])

points((d$cov_mean_diff)[ps], d$cov_cv_diff[ps],
       pch = 2,
       cex = cex_level, col = adjustcolor(collst[4], alpha_level))

polygon(c(xsq, rev(xsq)), c(pred_cngd[,c(3)], rev(pred_cngd[,c(4)])), col = adjustcolor(collst[4], 0.5), border = NA)
lines(xsq, pred_cngd[,1], lty = c(1), col = collst[4], lwd = 1.8)
#matlines(xsq, pred_cngd[,c(1,4,5)], lty = c(1,2,2), col = collst[4])

#legend("topright", c("Same Surveyor", "Different Surveyors"),
#       lty = 1, pch = c(1,2), col = collst[3:4], bty = "n")
title("C. Non-Graminoids", adj = 0)

mtext("Mean non-zero species-level cover", 1, line = 2.4)
mtext("CV of non-zero species-level cover", 2, line = 2.4, adj = -0.3)
dev.off()

