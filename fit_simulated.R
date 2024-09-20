require(BayesianTools)
rm(list=ls())

doplot = TRUE

source("~/Dropbox/Rfunctions/logit_funs.R")
source("~/Dropbox/Rfunctions/truncnorm_funs.R")

### make data
set.seed(1234)

## parameters
b0 = 0.05
b1 = 0.2 # observation error parameters
n = 100  # sample size

## hyperparamters
xT_mu = logit(0.2)
xT_sd = abs(logit(0.2))/2

## true states
xT = rnorm(n, xT_mu, xT_sd)
x_sd = b0 + b1*ilogit(xT) # observation error sd

## observed states
xO = cbind(rnorm(n, xT, x_sd),
           rnorm(n, xT, x_sd))
xmu = rowMeans(xO)

if(doplot) {
  #dev.new()
  hist(ilogit(xT), breaks = 20)
  plot(ilogit(xT), x_sd); abline(a=b0, b=b1)
  matplot(ilogit(xT), ilogit(xO))
}


### fit model
nparam = 2 # number of parameters
# param = c(b0, b1, xT_mu, xT_sd, xT)
param = c(b0, b1,  xT)

## likelihood function
likelihood <- function(param){
  # estimated parameters
  b0 = param[1]
  b1 = param[2]
  
  # estimated hyperparameters
  xT_mu = xT_mu#param[3]
  xT_sd = xT_sd#param[4]
  
  # estimated states
  xE = param[(nparam+1):length(param)]
  x_sd_est = b0+ilogit(xE)*b1
  
  LL1 = sum(dnorm(xO[,1], xE, x_sd_est, log = TRUE))
  LL2 = sum(dnorm(xO[,2], xE, x_sd_est, log = TRUE))
  
  llObservation =  LL1+LL2
  
  llRandomeffects = sum(dnorm(xE, xT_mu, xT_sd, log = TRUE))

  return(llObservation+llRandomeffects)
}

## set up Bayesian run
#c(b0, b1, xT_mu, xT_sd, xT)
setup <- createBayesianSetup(likelihood = likelihood,
                             lower = c(0, 0, rep(-4, n)),
                             upper = c(1, 2, rep(4, n)))

## run MCMC and save outputs
nmcmc = 1e5
settings <- list(iterations = nmcmc, consoleUpdates=1e1)
res <- runMCMC(bayesianSetup = setup, settings = settings)

smpout = getSample(res, start = (nmcmc/3)/2, parametersOnly=FALSE)
xE_mu = colMeans(smpout[,1:nparam])
xE_sd = apply(smpout[,1:nparam], 2, sd)
gelman_out = gelmanDiagnostics(res)
gelman_out

## plot diagnostics for a single MCMC run
plot(res, start = (nmcmc/3)/2, whichParameters = 1:nparam)
round(colMeans(smpout[,1:nparam]),3)
c(b0, b1, xT_mu, xT_sd)

hist(smpout[,1], xlim = c(0,1)); abline(v = b0, lty = 2, col = 2, lwd = 2)
hist(smpout[,2], xlim = c(0,1)); abline(v = b1, lty = 2, col = 2, lwd = 2)
#hist(smpout[,3], xlim = c(0,1)); abline(v = xT_mu, lty = 2, col = 2, lwd = 2)
#hist(smpout[,4], xlim = c(0,1)); abline(v = xT_sd, lty = 2, col = 2, lwd = 2)

xTest = colMeans(smpout[,(nparam+1):(nparam+n)])
xTest_sd = apply(smpout[,(nparam+1):(nparam+n)],2,function(x) quantile(x, pnorm(c(-1,1))))
plot(xT, xTest); abline(a = 0, b = 1, lty = 2)
segments(xT, xTest_sd[1,], xT, xTest_sd[2,])

# plot estimate for observation i
i = sample(n, 1)
hist(ilogit(smpout[,(nparam+1):(nparam+n)][,i]), breaks = 20, main = "xT[i]", xlim = c(0,4)); abline(v=xT[i], col = 2, lwd=2)
