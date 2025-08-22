require(BayesianTools)
rm(list=ls())

# make data
b0 = 0.01
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

moddat = data.frame(xtrue = xT, xO1 = xO[,1], xO2 = xO[,2], xOmu = rowMeans(xO))
xmu = moddat$xOmu
deltax = moddat$xO1-moddat$xO2
x_sd_est = sqrt((deltax)^2/2)

xmu_sd_est = x_sd_est/sqrt(2)

lvl = cut(xmu, 50) # cut(xT, 50) # 
xmu_est = tapply(xmu, lvl, mean)
xsd_st = sqrt(tapply(x_sd_est^2, lvl, mean))
xsd_st_strue = sqrt(tapply(x_sd^2, lvl, mean))
xmu_true = tapply(xT, lvl, mean)

plot(xmu_est, xsd_st)
points(xmu_true,xsd_st_strue, col=2)

plot(xmu_est, (xsd_st-b0-b1*xmu_est))


#abline(a=b0, b=b1)

cfs = coef(gls(xsd_st~xmu_est, weights = varExp(form = ~fitted(.)), na.action = na.exclude))
cfout[i,] = cfs

pars = c(0.1, 0.05)
xT = rnorm(n, 0.5, 0.1) #ilogit(rnorm(n, 0, 1))
x_sd = b0 + b1*xT
xO = cbind(rnorm(n, xT, x_sd),
             rnorm(n, xT, x_sd))
moddat = data.frame(xtrue = xT, xO1 = xO[,1], xO2 = xO[,2], xOmu = rowMeans(xO))


likelihood <- function(pars){
  b0 = pars[1]
  b1 = pars[2]
  xmu = moddat$xOmu
  
  deltax = moddat$xO1-moddat$xO2
  x_sd_est = sqrt((deltax)^2/2)
  x_est = (x_sd_est-b0)/b1
  
  
  # variance of xOmu = x_sd/sqrt(2)
  #x_sd_reg = b0 + b1*xmu
  #llObservation = dnorm(x_sd_reg, x_sd_est, x_sd_est)
  #lvl = cut(xmu, 20)
  #xmu_est = tapply(xmu, lvl, mean)
  #xsd_st = sqrt(tapply(x_sd_est^2, lvl, mean))
  
  #lm(xsd_st~xmu_est)
  llObservation = sum(dnorm(xmu, x_est, x_sd_est/sqrt(2)))
  #llObservation = sum(dnorm(moddat$xO1-xmu, sd = x_sd, log = TRUE))+
  #  sum(dnorm(moddat$xO2-xmu, sd = x_sd, log = TRUE))
  return(llObservation)
}


setup <- createBayesianSetup(likelihood = likelihood, lower = c(0,0), upper = c(1, 1))
settings <- list(iterations = 3e4,burnin = 2e4,consoleUpdates=1e3)
res <- runMCMC(bayesianSetup = setup, settings = settings)

smpout = getSample(res)
colMeans(smpout)
pars












#hist(cfout[,1], breaks = 20); abline(v=b0)
#hist(cfout[,2], breaks = 20); abline(v=b1)

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
