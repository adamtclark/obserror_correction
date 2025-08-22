require(BayesianTools)
rm(list=ls())

# make data
b0 = 0
b1 = 0.2

n = 1e3
pars = c(b0, b1)
niter = 1
cfout = matrix(nrow = niter, ncol = 2)

xT_mu = 0.5
xT_sd = 0.1

xT = rnorm(n, xT_mu, xT_sd)
hist(xT)
x_sd = b0 + b1*xT
plot(xT, x_sd)

xO = cbind(rnorm(n, xT, x_sd),
             rnorm(n, xT, x_sd))

matplot(xT, xO, pch = c(1,2))

xmu = rowMeans(xO)
plot(xT, xmu)

deltax = xO[,1]-xO[,2]
x_sd_est = sqrt((deltax)^2/2)
moddat = data.frame(xtrue = xT, xO1 = xO[,1], xO2 = xO[,2], xOmu = xmu,
                    xO_sd = x_sd_est)
#param = pars

require(lognorm)

likelihood <- function(param){
  #b0 = param[1]
  b1 = param[1]

  xmu = moddat$xOmu
  xO_sd = moddat$xO_sd
  
  x_sd_est = b0 + b1*xmu #moddat$xO_sd
  #x_sd_est = b0 + b1*(xT+rnorm(x_sd))
  #x_sd_est = b0 + b1*xT+b1*rnorm(x_sd)
  #x_sd_est = b0 + cov(xmu, s_sd)/(sigmaT^2+x_sd^2)*xmu
  
  #b1* = cov(xmu, x_sd)/(sigmaT^2+x_sd^2)
  #1/b1* = sigmaT^2/cov(xmu, s_sd)+x_sd^2/cov(xmu, s_sd)
  #1/b1* = 1/b1+x_sd^2/cov(xmu, s_sd)
  #b1* = 1/(1/b1+x_sd^2/cov(xmu, s_sd))
  
  # get var of x_sd
  lvl = cut(xmu, 20)
  xmu_sq = tapply(xmu, lvl, mean)
  xsd_var_sq = tapply(xO_sd^2, lvl, sd)
  #plot(xmu_sq, xsd_var_sq)
  mod = lm(xsd_var_sq~xmu_sq)
  #abline(mod)
  var_x0_var = predict(mod, newdata = data.frame(xmu_sq = xmu))^2
  #var_x0_var[var_x0_var<0] = min(var_x0_var[var_x0_var>0])
  #var(x^2) = 2*mu^4*sigma^2+2*sigma^4
  #var_x0_var = 2*x_sd_est^4+var_x0_sd+2*var_x0_sd^2
  
  #llObservation = sum(dnorm(xO_sd^2-x_sd_est^2, sd = sqrt(var_x0_var), log = TRUE))
  
  lmpar = getParmsLognormForMoments(mean = x_sd_est^2,
                                    var = var_x0_var)
  llObservation = sum(dlnorm(xO_sd^2, meanlog = lmpar[,1],
                             sd = sqrt(lmpar[,2]), log = TRUE))
  
  
  #llObservation = sum(dnorm(moddat$xO1-xmu, sd = x_sd_est, log = TRUE))+
  #  sum(dnorm(moddat$xO2-xmu, sd = x_sd, log = TRUE))
  return(llObservation)
}


#setup <- createBayesianSetup(likelihood = likelihood, lower = c(0,0), upper = c(1, 1))
setup <- createBayesianSetup(likelihood = likelihood, lower = c(0), upper = c(1))
settings <- list(iterations = 3e4,burnin = 1e4,consoleUpdates=1e2)
res <- runMCMC(bayesianSetup = setup, settings = settings)
plot(res, start = 100)

smpout = getSample(res, start = 100,parametersOnly=FALSE)
round(colMeans(smpout[,1:2]),3)
pars

hist(smpout[,1], breaks = 20); abline(v=b0, col = 2, lwd=2)
#hist(smpout[,2], breaks = 20); abline(v=b1, col = 2, lwd=2)
plot(smpout[,"par 1"], smpout[,"Lposterior"])
abline(v = b1, lty = 2)







b0^2+# check SD
sd_est0 = b0 + b1*moddat$xOmu

rtmp = rnorm(n, 0, x_sd)
sd_est = b0 + b1*(moddat$xtrue+rtmp/sqrt(2))

sd_est = b0 + b1*moddat$xtrue+b1*rtmp/sqrt(2)
sd_est = x_sd+b1*rtmp/sqrt(2)

x_sd_corrected = sd_est-b1*rtmp/sqrt(2)
x_sd_corrected = sd_est-b1*rtmp/sqrt(2)


plot(x_sd, x_sd_corrected); abline(a=0, b=1)

#x_sd = b0 + b1*xT
#sd_est = x_sd + b1*rnorm(n, 0, x_sd/sqrt(2))

plot(x_sd, sd_est); abline(a=0, b=1, lty=3)
points(x_sd, sd_est0, col = 2)


tmp = (xmu-xT)#/(x_sd/sqrt(2))
sq = seq(0,max(xT),length=20)
lvls=cut(xT, breaks = sq)
tmp2 = tapply(tmp, lvls, sd)
plot(sq[-1], tmp2)
lines(sq[-1], (b0+b1*sq[-1])/sqrt(2))


plot(xT, tmp)

coef(lm(x_sd~tmp))




lvls=cut(xT, breaks = sq)
tmp = cbind(tapply(x_sd, lvls, sd),
      tapply(sd_est0, lvls, sd))
matplot(sq[-1], tmp)





#sd_est-b1*rnorm(n, 0, x_sd/sqrt(2)) = b0 + b1*moddat$xtrue
#sd_est-b1*rnorm(n, 0, sd_est/sqrt(2)) = b0 + b1*moddat$xtrue










parssummary(res)
plot(res, burnin = 2e4)
marginalPlot(res, prior = T)
