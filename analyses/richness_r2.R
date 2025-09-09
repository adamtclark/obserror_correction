setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
#require(brms)

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.3
cex_level = 0.8
xsq = seq(0, 1, length = 1e3)

d = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
d = d[d$trt=="Control",]

###################
### number of species
hillfun = function(x, q = 0) {
  x[x==0] = NA
  if(is.null(dim(x))) {
    x = matrix(x)
  }
  if(nrow(x)>1) {
    p = x/rep(colSums(x, na.rm=TRUE), each = nrow(x))
    if(q == 0) {
      D = colSums(x>0, na.rm=TRUE)
    } else if(q == 1) {
      D = exp(-colSums(p*log(p), na.rm = TRUE))
    } else {
      D = colSums(p^q,na.rm=TRUE)^(1/(1-q))
    }
    D[!is.finite(D)] = NA
    D[apply(x, 2, function(y) all(is.na(y)))] = NA
    return(D)
  } else {
    return(NA)
  }
}

## richness
d_resurvey_ag = with(d_resurvey[!d_resurvey$species%in%c("bare ground", "litter", "cryptogam"),], aggregate(list(div_1 = cover_1,
                                                                                                                 div_2 = cover_2,
                                                                                                                 div_3 = cover_3,
                                                                                                                 div_4 = cover_4),
                                                                                                            by = list(SITE_CODE = SITE_CODE, trt = trt,
                                                                                                                      first_nutrient_year = first_nutrient_year,
                                                                                                                      block = block, plot = plot),
                                                                                                            FUN = function(x) hillfun(x, q = 0)))
agfun = function(q) {
  with(d_resurvey[!d_resurvey$species%in%c("bare ground", "litter", "cryptogam"),], aggregate(list(div_1 = cover_1,
             div_2 = cover_2,
             div_3 = cover_3,
             div_4 = cover_4),
        by = list(SITE_CODE = SITE_CODE, trt = trt,
                  first_nutrient_year = first_nutrient_year,
                  block = block, plot = plot),
        FUN = function(x) hillfun(x, q = q)))
}

d_resurvey_ag_RI = agfun(0)

collst = viridis(3)

pdf("figures/diversity_estimates.pdf", width = 10, height = 4)
par(mfcol = c(2,3), mar = c(4,4,2,2))
plot(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_2, xlab = "richness, first survey", ylab = "richness, repeat surveys", col = collst[1])
points(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_3, col = collst[2])
points(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_4, col = collst[3])
abline(a=0, b = 1, lty = 2)
slopes = plot_sites_fun(d_resurvey_ag_RI)
legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
abline(v = mean(slopes, na.rm=TRUE), lty = 2)


## shannon
d_resurvey_ag_SH = agfun(1)

collst = viridis(3)
plot(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_2, xlab = "Shannon, first survey", ylab = "Shannon, repeat surveys", col = collst[1])
points(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_3, col = collst[2])
points(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_4, col = collst[3])
abline(a=0, b = 1, lty = 2)
slopes = plot_sites_fun(d_resurvey_ag_SH)
legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
abline(v = mean(slopes, na.rm=TRUE), lty = 2)

## Simpson
d_resurvey_ag_SI = agfun(2)

collst = viridis(3)
plot(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_2, xlab = "Shannon, first survey", ylab = "Shannon, repeat surveys", col = collst[1])
points(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_3, col = collst[2])
points(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_4, col = collst[3])
abline(a=0, b = 1, lty = 2)
slopes = plot_sites_fun(d_resurvey_ag_SI)
legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
abline(v = mean(slopes, na.rm=TRUE), lty = 2)

# note - higher error, due to impact of cover?
dev.off()