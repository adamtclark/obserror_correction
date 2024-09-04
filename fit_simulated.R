require(BayesianTools)
rm(list=ls())

set.seed(240904)
doplot = FALSE

### make data
## parameters
b1 = 0.2 # observation error parameters
n = 100 # sample size

## hyperparamters
xT_mu = 0.2
xT_sd = 0.1

## true states
xT = rlnorm(n, xT_mu, xT_sd)
x_sd = b1*xT # observation error sd

## observed states
xO = cbind(rnorm(n, xT, x_sd),
           rnorm(n, xT, x_sd))
xmu = rowMeans(xO)
#x_sd_est = sqrt((xO[,1]-xO[,2])^2/2)

#plot(xmu, x_sd_est)
#require(lmodel2)
#tmp[i] = lmodel2(x_sd_est~xmu, range.y = "RMA", range.x = "RMA")$regression.results[4,3]


## standardised estimates
xO_std_est = xO/xmu

if(doplot) {
  plot(xT, x_sd, log = "xy")
  
  bins = cut(xT, 20)
  sqrt(mean(apply(xO_std_est, 1, var)))
  plot(tapply(xT, bins, mean), tapply(apply(xO_std_est, 1, var), bins, function(x) sqrt(mean(x))))
  
  matplot(xT, xO, pch = c(1,2), log = "xy")
  plot(xT, xmu)
  
  op=par(mfrow=c(2,2))
  mod = lm(log(xT) ~ log(xmu))
  summary(mod)
  plot(mod)
  par(op)
}


### fit model
nparam = 3 # number of parameters
# param = c(b1, xT_mu, xT_sd, xT)

## likelihood function
likelihood <- function(param){
  # estimated parameters
  x_sd_est = b1 = param[1]
  
  # estimated hyperparameters
  xT_mu = param[2]
  xT_sd = param[3]
  
  # estimated states
  xE = param[4:length(param)]
  
  llObservation =   sum(c(dnorm(xE-xO[,1], sd = xE*x_sd_est, log = TRUE), 
                          dnorm(xE-xO[,2], sd = xE*x_sd_est, log = TRUE)))
  llRandomeffects = sum(dlnorm(xE, xT_mu, xT_sd, log = TRUE))
  
  return(llObservation+llRandomeffects)
}

## set up Bayesian run
setup <- createBayesianSetup(likelihood = likelihood,
                             lower = c(0, 0, 0, rep(0, n)),
                             upper = c(2, 5, 2, rep(2, n)))

## repeat fitting niter times
niter = 1
datout = matrix(nrow = niter, ncol = nparam*3+3)
stateout = matrix(nrow = niter, ncol = n)
truestates = matrix(nrow = niter, ncol = n)

j = 1; i =1
for(j in 1:niter) {
  ## generate data and get "naive" estimates
  
  vary_states = FALSE
  if(vary_states) {
    # true states
    xT = rnorm(n, xT_mu, xT_sd)
    x_sd = b1*xT # observation error sd
    truestates[j,] = xT
    
    # observed states
    xO = cbind(rnorm(n, xT, x_sd),
               rnorm(n, xT, x_sd))
    xmu = rowMeans(xO)
    
    # standardised estimates
    xO_std_est = xO/xmu
    ps = 1:n
  } else {
    ps = sample(1:n, rep = TRUE)
  }
  
  ## run MCMC and save outputs
  settings <- list(iterations = 2e5, consoleUpdates=1e1)
  res <- runMCMC(bayesianSetup = setup, settings = settings)
  
  smpout = getSample(res, start = 6e4, parametersOnly=FALSE)
  xE_mu = colMeans(smpout[,1:nparam])
  xE_sd = apply(smpout[,1:nparam], 2, sd)
  
  # diagnostics
  datout[j,1:nparam] = gelmanDiagnostics(res, start = 5e4, whichParameters = 1:nparam)$psrf[,1]
  
  # parameters
  datout[j,1:nparam+nparam] = xE_mu
  datout[j,1:nparam+nparam*2] = xE_sd
  
  # naive estimates
  parsest = sqrt(mean(apply(xO_std_est[ps,], 1, var)))
  datout[j,nparam*3+1] = parsest # naive parameter estimate, b1
  datout[j,nparam*3+2] = mean(xmu[ps]) # naive parameter estimate, xT_mu
  datout[j,nparam*3+3] = sd(xmu[ps]) # naive parameter estimate, xT_sd
  
  # states
  stateout[j,] = colMeans(smpout[,(nparam+1):(nparam+n)])
  
  print(j/niter)
}

# save(list = c("datout", "stateout", "truestates"), file = "output/mcmcout_fixed.rda")
# load("output/mcmcout.rda")
# load("output/mcmcout_fixed.rda")

## check hyperparameter likelihood
if(FALSE) {
  whichrun = sample(dim(smpout)[1], 1)
  xE = smpout[,(nparam+1):(nparam+n)][whichrun,]
  xmusq = seq(0, 1, length = 1000)
  tmp = numeric(1000)
  for(i in 1:1000) {
    tmp[i] = likelihood(c(b1, xmusq[i], xT_sd, xE))
      #sum(dnorm(xE-xmusq[i], 0.1, log = TRUE))
  }
  plot(xmusq, tmp, xlab = "xT_mu", ylab = "LL")
  abline(v=xT_mu, lty=2)
}

## plot diagnostics for a single MCMC run
plot(res, start = 6e4, whichParameters = 1:nparam)
round(colMeans(smpout[,1:nparam]),3)
c(b1, xT_mu, xT_sd)

hist(smpout[,1], xlim = c(0,1)); abline(v = b1, lty = 2, col = 2, lwd = 2)
hist(smpout[,2], xlim = c(0,1)); abline(v = xT_mu, lty = 2, col = 2, lwd = 2)
hist(smpout[,3], xlim = c(0,1)); abline(v = xT_sd, lty = 2, col = 2, lwd = 2)

xTest = colMeans(smpout[,(nparam+1):(nparam+n)])
xTest_sd = apply(smpout[,(nparam+1):(nparam+n)],2,function(x) quantile(x, pnorm(c(-1,1))))
plot(xT, xTest); abline(a = 0, b = 1, lty = 2)
segments(xT, xTest_sd[1,], xT, xTest_sd[2,])
gelmanDiagnostics(res, start = 100, whichParameters = 1:nparam)$psrf[,1]

# plot estimate for observation i
i = sample(n, 1)
hist(smpout[,(nparam+1):(nparam+n)][,i], breaks = 20, main = "xT[i]", xlim = c(0,1)); abline(v=xT[i], col = 2, lwd=2)

if(FALSE) {
  ### check for bias
  ## parameters
  par(mar=c(4,4,2,2), mfrow=c(2,1))
  hist(datout[,1+nparam], breaks = 20, main = "b1", xlim = c(b1/2,b1*2)); abline(v=b1, col = 2, lwd=2)
  hist(datout[,1+3*nparam], breaks = 20, main = "b1_naive", xlim = c(b1/2,b1*2)); abline(v=b1, col = 2, lwd=2)
  
  ## hyperparameters
  hist(datout[,2+nparam], breaks = 20, main = "xT_mu", xlim = c(0,1)); abline(v=xT_mu, col = 2, lwd=2)
  hist(datout[,2+3*nparam], breaks = 20, main = "xT_mu_naive", xlim = c(0,1)); abline(v=xT_mu, col = 2, lwd=2)
  
  hist(datout[,3+nparam], breaks = 20, main = "xT_sd", xlim = c(0,1)); abline(v=xT_sd, col = 2, lwd=2)
  hist(datout[,3+3*nparam], breaks = 20, main = "xT_sd_naive", xlim = c(0,1)); abline(v=xT_sd, col = 2, lwd=2)
  
  ## diagnostics
  par(mar=c(4,4,2,2), mfrow=c(1,1))
  hist(datout[,1], breaks = 20, main = "Gelman b1"); abline(v=1, col = 2, lwd=2)
  hist(datout[,2], breaks = 20, main = "Gelman xT_mu"); abline(v=1, col = 2, lwd=2)
  hist(datout[,2], breaks = 20, main = "Gelman xT_sd"); abline(v=1, col = 2, lwd=2)
  
  
  ## states
  if(vary_states) {
    plot(truestates, stateout); abline(a = 0, b = 1, lty = 2, col = 2, lwd = 2)
  } else {
    xTest = colMeans(stateout)
    xTest_sd = sqrt(colMeans((stateout*datout[,1+nparam])^2)) # apply(stateout,2,function(x) quantile(x, pnorm(c(-1,1))))
    
    plot(xT, xTest); abline(a = 0, b = 1, lty = 2)
    segments(xT, xTest-xTest_sd, xT, xTest+xTest_sd)
    
    points(xT, xmu, col = 2)
    
    deltax = xO[,1]-xO[,2]
    x_sd_est_raw = sqrt((deltax)^2/2)
    
    segments(xT, xmu+x_sd_est_raw, xT, xmu-x_sd_est_raw, col = 2)
  }
  
  # plot individual states
  i = 41
  hist(stateout[,i], breaks = 20, main = "xT[i]", xlim = c(0,1)); abline(v=xT[i], col = 2, lwd=2)
}
