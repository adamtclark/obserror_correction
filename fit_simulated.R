require(BayesianTools)
rm(list=ls())

# make data
b0 = 0
b1 = 0.2

n = 50
pars = c(b0, b1)
niter = 1
cfout = matrix(nrow = niter, ncol = 2)

xT_mu = 0.5
xT_sd = 0.1

xT = rnorm(n, xT_mu, xT_sd)
#hist(xT)
x_sd = b0 + b1*xT
plot(xT, x_sd, log = "xy")

xO = cbind(rnorm(n, xT, x_sd),
             rnorm(n, xT, x_sd))
xmu = rowMeans(xO)
xsigma = sqrt((xO[,1]-xO[,2])^2/2)

bins = cut(xT, 20)
xO_std = xO
xO_std = xO_std/xmu

sqrt(mean(apply(xO_std, 1, var)))

plot(tapply(xT, bins, mean), tapply(apply(xO_std, 1, var), bins, function(x) sqrt(mean(x))))


matplot(xT, xO, pch = c(1,2), log = "xy")
#plot(xT, xmu)
#plot(lm(log(xT) ~ log(xmu)))


deltax = xO[,1]-xO[,2]
x_sd_est = sqrt((deltax)^2/2)
moddat = data.frame(xtrue = xT, xO1 = xO[,1], xO2 = xO[,2], xOmu = xmu,
                    xO_sd = x_sd_est)

# fit model
# param = c(b0, b1, xT_mu, xT_sd, rnorm(n, xT_mu, xT_sd))
likelihood <- function(param){
  b0 = param[1]
  b1 = param[2]

  xT_mu = param[3]
  xT_sd = param[4]
  
  xE = param[5:length(param)]
  
  x_sd_est = b0 + b1*xE
  
  llObservation =   sum(c(dnorm(xE-moddat$xO1, sd = x_sd_est, log = TRUE), 
                        dnorm(xE-moddat$xO2, sd = x_sd_est, log = TRUE)))
  #llObservation =   sum(dnorm(xE-moddat$xOmu, sd = x_sd_est/sqrt(2), log = TRUE))

  llRandomeffects = sum(dnorm(xE-xT_mu, xT_sd, log = TRUE))

  return(llObservation+llRandomeffects)
}


#setup <- createBayesianSetup(likelihood = likelihood, lower = c(0,0), upper = c(1, 1))
setup <- createBayesianSetup(likelihood = likelihood,
                             lower = c(0, 0, 0, 0, rep(0, n)),
                             upper = c(1, 1, 1, 1, rep(1, n)))
nchains = 4

niter = 100
datout = matrix(nrow = niter, ncol = 12)
for(j in 1:niter) {
  sV = NULL
  for(i in 1:nchains) {
    ps = sample(1:length(xmu), rep = TRUE)
    parsest = pmax(coef(lm(xsigma[ps]~xmu[ps])), 1e-3)
    
    sV = rbind(sV, c(
      parsest, # b0, b1
      mean(xmu[ps]), #xT_mu
      mean(xsigma[ps]), #xT_sd
      xmu+rnorm(length(xmu), 0, (parsest[1]+parsest[2]*xmu)/sqrt(2))
    ))
  }
  settings <- list(iterations = 5e4, burnin = 3e4, consoleUpdates=1e1, startValue=sV)
  
  res <- runMCMC(bayesianSetup = setup, settings = settings)
  
  smpout = getSample(res, start = 100, parametersOnly=FALSE)
  xE_mu = colMeans(smpout[,1:4])
  xE_sd = apply(smpout[,1:4], 2, sd)
  
  datout[j,1:4] = gelmanDiagnostics(res, start = 100, whichParameters = 1:4)$psrf[,1]
  datout[j,5:8] = xE_mu[1:4]
  datout[j,9:12] = xE_sd[1:4]
  
  print(j/niter)
}

# check for bias
par(mar=c(4,4,2,2))
hist(datout[,5], breaks = 20, main = "b0"); abline(v=b0, col = 2, lwd=2)
hist(datout[,6], breaks = 20, main = "b1"); abline(v=b1, col = 2, lwd=2)
hist(datout[,7], breaks = 20, main = "xT_mu"); abline(v=xT_mu, col = 2, lwd=2)
hist(datout[,8], breaks = 20, main = "xT_sd"); abline(v=xT_sd, col = 2, lwd=2)

hist(datout[,1], breaks = 20, main = "Gelman b0"); abline(v=1, col = 2, lwd=2)
hist(datout[,2], breaks = 20, main = "Gelman b1"); abline(v=1, col = 2, lwd=2)
hist(datout[,3], breaks = 20, main = "Gelman xT_mu"); abline(v=1, col = 2, lwd=2)
hist(datout[,4], breaks = 20, main = "Gelman xT_sd"); abline(v=1, col = 2, lwd=2)







if(FALSE) {
  plot(res, start = 100, whichParameters = 1:2)
  plot(res, start = 100, whichParameters = 3:4)
  plot(res, start = 100, whichParameters = 5:6)
  
  round(colMeans(smpout[,1:2]),3)
  pars
  
  hist(smpout[,1], xlim = c(0,1)); abline(v = b0, lty = 2, col = 2, lwd = 2)
  hist(smpout[,2], xlim = c(0,1)); abline(v = b1, lty = 2, col = 2, lwd = 2)
  
  hist(smpout[,3], xlim = c(0,1)); abline(v = xT_mu, lty = 2, col = 2, lwd = 2)
  hist(smpout[,4], xlim = c(0,1)); abline(v = xT_sd, lty = 2, col = 2, lwd = 2)
  
  xTest = colMeans(smpout[,5:(5+n-1)])
  xTest_sd = apply(smpout[,5:(5+n-1)],2,sd)
  plot(xT, xTest); abline(a = 0, b = 1, lty = 2)
  segments(xT, xTest+xTest_sd, xT, xTest-xTest_sd)
  
  #hist(smpout[,1], breaks = 20); abline(v=b0, col = 2, lwd=2)
  #hist(smpout[,2], breaks = 20); abline(v=b1, col = 2, lwd=2)
  plot(smpout[,"par 1"], smpout[,"Lposterior"])
  abline(v = b1, lty = 2)
}