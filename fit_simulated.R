require(BayesianTools)
rm(list=ls())

set.seed(240905)
doplot = TRUE

### truncated normal functions
#x = c(0.2, 0.6, 0.9); mu = 0.3; sigma = 0.5; a = 0; b = 1; type = c("density", "mean")
truncnorm = function(x = NULL, mu = 0.5, sigma = 1, a = 0, b = 1, type = c("density")) {
  #a and b are min/max values
  #mu and sigma the mean and std. dev.
  #x is value for which to calculate the function
  #type can be density, cumdensity, mean, and/or variance
  
  xi = (x-mu)/sigma
  alpha = (a-mu)/sigma
  beta = (b-mu)/sigma
  Z = pnorm(beta) - pnorm(alpha)
  
  out = NULL
  n = 1
  if(sum(type=="density")>0) {
    out[[n]] = dnorm(xi)/(sigma*Z)
    names(out)[n] = "density"
    n = n+1
  }
  if(sum(type=="cumdensity")>0) {
    out[[n]] = (pnorm(xi)-pnorm(alpha))/Z
    names(out)[n] = "cumdensity"
    n = n+1
  }
  if(sum(type=="mean")>0) {
    out[[n]] = mu + sigma*(dnorm(alpha)-dnorm(beta))/Z
    names(out)[n] = "mean"
    n = n+1
  }
  if(sum(type=="variance")>0) {
    out[[n]] = sigma^2*(1-(beta*dnorm(beta)-alpha*dnorm(alpha))/Z-((dnorm(alpha)-dnorm(beta))/Z)^2)
    names(out)[n] = "variance"
    n = n+1
  }
  
  return(out)
}

# check
if(FALSE) {
  # density
  dx = 1e-3
  xseq = seq(0, 1, by = dx)
  out_dens = truncnorm(x = xseq, mu = 0.2, sigma = 0.5, a = 0, b = 1,
                  type = "density")$density
  plot(xseq, out_dens, type = "l")
  sum(out_dens*dx)
  
  # cum density
  out_cdens = truncnorm(x = xseq, mu = 0.2, sigma = 0.5, a = 0, b = 1,
                  type = "cumdensity")$cumdensity
  plot(xseq, out_cdens, type = "l")
  plot(out_cdens, cumsum(out_dens*dx)); abline(a=0, b=1, lty=2, col = "blue")
  
  # mean
  out_mean = truncnorm(x = xseq, mu = 0.2, sigma = 0.5, a = 0, b = 1,
                        type = "mean")$mean
  c(out_mean, sum(xseq*(out_dens/sum(out_dens))))
  
  # variance
  out_var = truncnorm(x = xseq, mu = 0.2, sigma = 0.5, a = 0, b = 1,
                       type = "variance")$variance
  c(out_var, sum((xseq-out_mean)^2*(out_dens/sum(out_dens))))
}

## make inverse functions?
## make random number sampler


### make data
## parameters
b1 = 0.2 # observation error parameters
n = 100  # sample size

## hyperparamters
xT_mu = 0.01
xT_sd = 0.15

## true states
xT = rlnorm(n, xT_mu, xT_sd)
x_sd = b1*xT # observation error sd

## observed states
xO = cbind(rnorm(n, xT, x_sd),
           rnorm(n, xT, x_sd))
xmu = rowMeans(xO)

if(doplot) {
  hist(xT)
  plot(xT, x_sd)
  matplot(xT, xO)
}


### fit model
nparam = 3 # number of parameters
# param = c(b1, xT_mu, xT_sd, xT)

## likelihood function
likelihood <- function(param){
  # estimated parameters
  b1 = param[1]
  
  # estimated hyperparameters
  xT_mu = param[2]
  xT_sd = param[3]
  
  # estimated states
  xE = param[4:length(param)]
  x_sd_est = xE*b1
  
  llObservation =   sum(c(dnorm(xE-xO[,1], sd = x_sd_est, log = TRUE), 
                          dnorm(xE-xO[,2], sd = x_sd_est, log = TRUE)))
  llRandomeffects = sum(dlnorm(xE, xT_mu, xT_sd, log = TRUE))
  
  return(llObservation+llRandomeffects)
}

## set up Bayesian run
setup <- createBayesianSetup(likelihood = likelihood,
                             lower = c(0, 0, 0, rep(0, n)),
                             upper = c(2, 5, 2, rep(2, n)))

## run MCMC and save outputs
nmcmc = 2e5
settings <- list(iterations = nmcmc, consoleUpdates=1e1)
res <- runMCMC(bayesianSetup = setup, settings = settings)

smpout = getSample(res, start = (nmcmc/3)/2, parametersOnly=FALSE)
xE_mu = colMeans(smpout[,1:nparam])
xE_sd = apply(smpout[,1:nparam], 2, sd)
gelman_out = gelmanDiagnostics(res)

## plot diagnostics for a single MCMC run
plot(res, start = (nmcmc/3)/2, whichParameters = 1:nparam)
round(colMeans(smpout[,1:nparam]),3)
c(b1, xT_mu, xT_sd)

hist(smpout[,1], xlim = c(0,1)); abline(v = b1, lty = 2, col = 2, lwd = 2)
hist(smpout[,2], xlim = c(0,1)); abline(v = xT_mu, lty = 2, col = 2, lwd = 2)
hist(smpout[,3], xlim = c(0,1)); abline(v = xT_sd, lty = 2, col = 2, lwd = 2)

xTest = colMeans(smpout[,(nparam+1):(nparam+n)])
xTest_sd = apply(smpout[,(nparam+1):(nparam+n)],2,function(x) quantile(x, pnorm(c(-1,1))))
plot(xT, xTest); abline(a = 0, b = 1, lty = 2)
segments(xT, xTest_sd[1,], xT, xTest_sd[2,])

# plot estimate for observation i
i = sample(n, 1)
hist(smpout[,(nparam+1):(nparam+n)][,i], breaks = 20, main = "xT[i]", xlim = c(0,4)); abline(v=xT[i], col = 2, lwd=2)
