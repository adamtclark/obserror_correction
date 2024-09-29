require(BayesianTools)
#require(truncnorm)
#require(fitdistrplus)
rm(list=ls())
  
doplot = TRUE

source("~/Dropbox/Rfunctions/logit_funs.R")
#source("~/Dropbox/Rfunctions/truncnorm_funs.R")

### make data
set.seed(1234)

## hyperparamters
xT_mu = -2
xT_sd = 1

## parameters
b0 = -0.5
b1 = 0.5 # observation error parameters
n = 100  # sample size

## true states
xT = rnorm(n, mean = xT_mu, sd = xT_sd)
pT =  exp(xT)/(1 + exp(xT))
x_sd = exp(b0 + b1*xT)

## observed states
xO = cbind(rnorm(n, xT, x_sd),
           rnorm(n, xT, x_sd))
pO = exp(xO)/(1 + exp(xO))
xmu = rowMeans(xO)


if(doplot) {
  #dev.new()
  par(mfcol=c(3,2), mar=c(4,4,2,2))
  hist(pT, breaks = 20)
  plot(pT, x_sd)# abline(a=0, b=b1)
  matplot(pT, pO)
  
  hist(xT, breaks = 20)
  plot(xT, x_sd)# abline(a=0, b=b1)
  matplot(xT, xO)
}
#par(mfrow=c(1,1))
#dxest = (xT-xmu)/exp(b0+b1*xmu)
#plot(xT, dxest)
#plot(xT, dxest*exp(b0+b1*xmu)+xmu)

#xmu*exp(b0)*exp(b1*xT) = dxest

### Jesus fucking christ on a cracker -
# need to get the inverse of this for the sampler

### fit model
nparam = 2 # number of parameters
#param = c(b0, b1, xT_mu, xT_sd, xT)
#dxest = (xmu-xT)/xmu
#xT = xmu*(1-dxest)
#param = c(b0, b1, dxest)
param = c(b0, b1, xT)

#  xE = xmu/(1+1/sqrt(2)*dxest)
# 0 < xE < 10), upper = c(1,2)) 



## likelihood function
likelihood <- function(param){
  # estimated parameters
  b0 = param[1]
  b1 = param[2]
  
  # estimated hyperparameters
  xT_mu = mean(xmu)#xT_mu#param[3]
  xT_sd = sd(xmu)#xT_sd#param[4]
  
  # estimated states
  #dxest = param[-c(1:nparam)]
  #xE = xmu*(1-dxest)
  xE = (param[-c(1:nparam)])
  
  x_sd_est = exp(b0+xE*b1)
  
  LL1 = sum(dnorm(xO[,1], xE, x_sd_est, log = TRUE))
  LL2 = sum(dnorm(xO[,2], xE, x_sd_est, log = TRUE))
  
  llObservation =  LL1+LL2
  llRandomeffects = sum(dnorm(xE, xT_mu, xT_sd, log = TRUE))

  return(llObservation+llRandomeffects)
}


# try out likelihoods
if(doplot) {
  par(mfrow=c(1,1), mar=c(4,4,2,2))
  b0lst = seq(-1, 1, length=1e3)
  out = numeric(100)
  for(i in 1:length(b0lst)) {
    out[i] = likelihood(param = c(b0lst[i], b1, xT))
  }
  plot(b0lst, out, type = "l",
       ylim=quantile(out, c(0.5, 1))); abline(v = b0, lty =2)
  abline(v = b0lst[which.max(out)], col = 2, lty = 2)
  
  
  b1lst = seq(-1, 1, length=1e3)
  out = numeric(100)
  for(i in 1:length(b1lst)) {
    out[i] = likelihood(param = c(b0, b1lst[i], xT))
  }
  plot(b1lst, out, type = "l",
       ylim=quantile(out, c(0.5, 1))); abline(v = b1, lty =2)
  abline(v = b1lst[which.max(out)], col = 2, lty = 2)
}

## set up Bayesian run
#c(b0, b1, xT_mu, xT_sd, xT)
setup <- createBayesianSetup(likelihood = likelihood,
                             lower = c(-1, -1, rep(-5, n)),
                             upper = c(1, 1, rep(5, n)))
param = c(runif(2, -1, 1), runif(n, 0, 1))
likelihood(param)

## run MCMC and save outputs
nmcmc = 5e4
settings <- list(iterations = nmcmc, consoleUpdates=1e1)
res <- runMCMC(bayesianSetup = setup, settings = settings)


smpout = getSample(res, start = (nmcmc/3)/2, parametersOnly=FALSE)
xE_mu = colMeans(smpout[,1:nparam,drop = FALSE])
xE_sd = apply(smpout[,1:nparam,drop = FALSE], 2, sd)
gelman_out = gelmanDiagnostics(res)
gelman_out

## plot diagnostics for a single MCMC run
plot(res, start = (nmcmc/3)/2, whichParameters = 1:nparam)
round(colMeans(smpout[,1:nparam,drop = FALSE]),3)
c(b0, b1, xT_mu, xT_sd)

op = par(mfrow=c(2,1), mar=c(4,4,2,2))
hist(smpout[,1], xlim = c(-1,1), breaks = 20); abline(v = b0, lty = 2, col = 2, lwd = 2)
hist(smpout[,2], xlim = c(-1,1), breaks = 20); abline(v = b1, lty = 2, col = 2, lwd = 2)

#hist(smpout[,2], xlim = c(0,1), breaks = 20); abline(v = b1, lty = 2, col = 2, lwd = 2)
#hist(smpout[,3], xlim = c(0,1)); abline(v = xT_mu, lty = 2, col = 2, lwd = 2)
#hist(smpout[,4], xlim = c(0,1)); abline(v = xT_sd, lty = 2, col = 2, lwd = 2)

xTest = rowMeans(xmu*t(1-(smpout[,(nparam+1):length(param),drop=FALSE])))
xTest_sd = apply(xmu*t(1-(smpout[,(nparam+1):length(param),drop=FALSE])), 1, function(x) quantile(x, pnorm(c(-1,1))))

#del = rowMeans(t(1+1/sqrt(2)*(smpout[,(nparam+1):length(param),drop=FALSE])))
cor(xTest, xT)
cor(xmu, xT)
#xTest = rowMeans(xmu*(t(smpout[,(nparam+1):(nparam+n)])+1))
#xTest_sd = apply(xmu*(t(smpout[,(nparam+1):(nparam+n)])+1),1,function(x) quantile(x, pnorm(c(-1,1))))
plot(xT, xTest, ylim = range(c(xTest_sd[1,], xTest_sd[2,]))); abline(a = 0, b = 1, lty = 2)
segments(xT, xTest_sd[1,], xT, xTest_sd[2,])
par(op)

# plot estimate for observation i
par(mfrow=c(1,1), mar=c(4,4,2,2))
i = sample(n, 1)
hist(xmu[i]*(1-smpout[,(nparam+1):(nparam+n)][,i]), breaks = 20, main = "xT[i]", xlim = c(-5,5)); abline(v=xT[i], col = 2, lwd=2)




#### BRMS

require(brms)
xOdat = data.frame(plot = 1:n, obs = c(xO[,1], xO[,2]),
                   mu = rowMeans(xO), xT = xT)
bmout = brm(bf(obs ~ mu*(1 + obserr),
               obserr ~ 1|plot,
               b0+b1~1,
               nlf(sigma ~ b0+b1*(mu*(1 + obserr))),
               #sigma ~ xT,
               nl = TRUE),
            data = xOdat, family = gaussian())

# PROBLEM: get double predictions for sigma (one for each observation...)
# PROBLEM 2: can't get the obserror terms to be the right size...
# PROBLEM 3: distribution?...

bmout
plot(xT, (1+ranef(bmout)$plot[,1,])*xmu); abline(a=0, b=1, lty = 2)
plot(xT, xmu); abline(a=0, b=1, lty = 2)

plot(ranef(bmout)$plot[,1,], xT/xmu-1); abline(a=0, b=1, lty = 2)

sigout =  fitted(bmout, dpar = "sigma")[,1]
sigma_est = exp(fixef(bmout)["b0_Intercept",1]+fixef(bmout)["b1_Intercept",1]*xT)
plot(sigout, rep(sigma_est,2)); abline(a=0, b=1)
plot(sigout, rep(x_sd,2)); abline(a=0, b=1)

## how is sigma calcualted???


fixef(bmout)
#coef(bmout)











n = 100
p = runif(n)
x = log(p/(1-p))
p =  exp(x)/(1 + exp(x))

b1 = 0.24
sigma = b1*p
plot(x, sigma)
plot(p, sigma)

##############
sigma1 = b1*p
sigma2 = b1*exp(xT)/(1 + exp(xT))
plot(sigma1, sigma2); abline(a=0, b=1, lty=2)

xT = x
set.seed(123)
xobs1 = rnorm(n, xT, sigma2)
set.seed(123)
xobs2 = xT + b1*rnorm(n, sd = 1)*exp(xT)/(1 + exp(xT))
plot(xobs1, xobs2); abline(a=0, b=1, lty=2)

dx = xobs2-xT
set.seed(123); plot(dx,
                    b1*rnorm(n, sd = 1)*exp(xobs2-dx)/(1 + exp(xobs2-dx))); abline(a=0,b=1,lty=2)

set.seed(123); plot(xT,
                    xobs2-b1*rnorm(n, sd = 1)*exp(xobs2-dx)/(1 + exp(xobs2-dx))); abline(a=0,b=1,lty=2)


