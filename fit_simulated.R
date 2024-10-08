require(BayesianTools)
require(truncnorm)
require(fitdistrplus)
rm(list=ls())
  
doplot = TRUE

#source("~/Dropbox/Rfunctions/logit_funs.R")
#source("~/Dropbox/Rfunctions/truncnorm_funs.R")


### next steps
# try brms implementation
# try to add in probability of zeros

### notes?
# Forecast richness at a scale, and then estimate expected error effects and true changes. 
# Maybe add in prior from Emma? 
# Assumption: SAR, plus that abundance distribution is scale invariant 
# Add in distance matrix option (taxonomic, phylo, etc.)
# Use abundance distribution to estimate gains and losses


### make data
set.seed(1234)

## parameters
b0 = 0.01
b1 = 0.2 # observation error parameters
n = 50  # sample size

## hyperparamters
xT_mu = 0.2
xT_sd = 0.2

## true states
xT = rtruncnorm(n, a=0, b=1, mean = xT_mu, sd = xT_sd)
x_sd = b0 + b1*xT # observation error sd

## observed states
xO = cbind(rtruncnorm(n, a=0, b=1, xT, x_sd),
           rtruncnorm(n, a=0, b=1, xT, x_sd))
xmu = rowMeans(xO)

fitout = fitdist(xmu, "truncnorm", fix.arg=list(a=0, b=1),
        start = list(mean = mean(xmu), sd = sd(xmu)))
fitout

#need to think about this...


c(etruncnorm(a=0, b=1, xT_mu, xT_sd), mean(xmu))

if(doplot) {
  #dev.new()
  par(mfrow=c(3,1), mar=c(4,4,2,2))
  hist(xT, breaks = 20)
  plot(xT, x_sd); abline(a=b0, b=b1)
  matplot(xT, xO)
}


### fit model
nparam = 4 # number of parameters
#param = c(b0, b1, xT_mu, xT_sd, xT)
dxest = xT/xmu-1
param = c(b0, b1, xT_mu, xT_sd, dxest)

## likelihood function
likelihood <- function(param){
  # estimated parameters
  b0 = param[1]
  b1 = param[2]
  
  # estimated hyperparameters
  xT_mu = param[3]#xT_mu#mean(xmu)#
  xT_sd = param[4]#xT_sd#sd(xmu)#
  
  # estimated states
  xE = xmu*(1+param[(nparam+1):length(param)])
  x_sd_est = b0+xE*b1
  
  LL1 = sum(log(dtruncnorm(xO[,1], a=0, b=1, xE, x_sd_est)))
  LL2 = sum(log(dtruncnorm(xO[,2], a=0, b=1, xE, x_sd_est)))
  
  llObservation =  LL1+LL2
  llRandomeffects = sum(log(dtruncnorm(xE, a=0, b=1, xT_mu, xT_sd)))

  return(llObservation+llRandomeffects)
}

# try out likelihoods
if(doplot) {
  par(mfrow=c(2,1))
  b0lst = seq(0, 1, length=1e3)
  out = numeric(100)
  for(i in 1:length(b0lst)) {
    out[i] = likelihood(param = c(b0lst[i], b1, xT_mu, xT_sd, dxest))
  }
  plot(b0lst, out, ylim=quantile(out, c(0.5, 1)), type = "l"); abline(v = b0, lty =2)
  abline(v = b0lst[which.max(out)], col = 2, lty = 2)
  
  b1lst = seq(0, 1, length=1e3)
  out = numeric(100)
  for(i in 1:length(b1lst)) {
    out[i] = likelihood(param = c(b0, b1lst[i], xT_mu, xT_sd, dxest))
  }
  plot(b1lst, out, ylim=quantile(out, c(0.5, 1)), type = "l"); abline(v = b1, lty =2)
  abline(v = b1lst[which.max(out)], col = 2, lty = 2)
}

## set up Bayesian run
#c(b0, b1, xT_mu, xT_sd, xT)
setup <- createBayesianSetup(likelihood = likelihood,
                             #lower = c(0, 0), upper = c(1,2))
                             lower = c(0, 0, 0, 0, rep(-1, n)),
                             upper = c(1, 1, 1, 1, rep(1, n)))

## run MCMC and save outputs
nmcmc = 5e4
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

par(mfrow=c(4,1), mar=c(4,4,2,2))
hist(smpout[,1], xlim = c(0,0.1), breaks = 20); abline(v = b0, lty = 2, col = 2, lwd = 2)
hist(smpout[,2], xlim = c(0,1), breaks = 20); abline(v = b1, lty = 2, col = 2, lwd = 2)
hist(smpout[,3], xlim = c(0,1)); abline(v = xT_mu, lty = 2, col = 2, lwd = 2)
hist(smpout[,4], xlim = c(0,1)); abline(v = xT_sd, lty = 2, col = 2, lwd = 2)

par(mfrow=c(1,1))
xTest = rowMeans(xmu*(t(smpout[,(nparam+1):(nparam+n)])+1))
xTest_sd = apply(xmu*(t(smpout[,(nparam+1):(nparam+n)])+1),1,function(x) quantile(x, pnorm(c(-1,1))))
plot(xT, xTest, ylim = range(c(xTest_sd[1,], xTest_sd[2,]))); abline(a = 0, b = 1, lty = 2)
segments(xT, xTest_sd[1,], xT, xTest_sd[2,])


# plot estimate for observation i
i = sample(n, 1)
hist(xmu[i]*(1+smpout[,(nparam+1):(nparam+n)][,i]), breaks = 20, main = "xT[i]", xlim = c(0,1)); abline(v=xT[i], col = 2, lwd=2)

