library(dplyr)
library(rjags)
library(ecoforecastR)

datadir <- "/Users/rzabramoff/Documents/Forecasting/EF_tree_growth/"
precip <- read.csv(paste0(datadir,"buckhorn_daily_precip.csv"))
temps <- read.csv(paste0(datadir,"buckhorn_daily_temps.csv"))
growth <- read.csv(paste0(datadir,"buckhorn_growth_2007_2018.csv"))

RandomWalk_Indivs <- "
model{
x[1] ~ dnorm(x_ic,tau_ic)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)

for(t in 2:NT){ ## t = year
x[t]~dnorm(x[t-1],tau_add)  ## process model
}

for(i in 1:n){ ## i = individual trees
y[i]  ~ dnorm(x[time[i]],tau_obs)		        ## data model
}
}
"

y <- growth$HT_cm
time <- rep(1:10, 970)

data <- list(n=length(y),time=time,y=y,NT=10,x_ic=5,tau_ic=100,a_obs=1,r_obs=1,a_add=1,r_add=1)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)),tau_obs=5/var(y.samp))
}

j.model   <- jags.model(file = textConnection(RandomWalk_Indivs),
                        data = data,
                        inits = init,
                        n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_obs","tau_add"),
                            n.iter = 5000)

##evaluate model
plot(jags.out)
gelman.plot(jags.out)

##remove burnin
burnin = 500                                ## determine convergence
jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
gelman.plot(jags.burn)

acfplot(jags.burn)
effectiveSize(jags.burn)
summary(jags.burn)
out <- as.matrix(jags.burn)

## credible and prediction intervals
nsamp <- 5000
samp <- sample.int(nrow(out),nsamp)
xpred <- seq(1,10,length=10)  					## sequence of x values we're going to
npred <- length(xpred)				##      make predictions for
ypred <- matrix(0.0,nrow=nsamp,ncol=npred)	## storage for predictive interval
ycred <- matrix(0.0,nrow=nsamp,ncol=npred)	## storage for credible interval

for(g in seq_len(nsamp)){
  theta = out[samp[g],]
  ycred[g,] <- theta[grep("x",names(theta))]
  ypred[g,] <- rnorm(npred, ycred[g,], 1/sqrt(theta["tau_obs"]))
}

ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))  ## credible interval and median
pi <- apply(ypred,2,quantile,c(0.025,0.975))		## prediction interval

plot(data$time,data$y,cex=0.5)
lines(xpred,ci[1,],col=3,lty=2)	## lower CI
lines(xpred,ci[2,],col=3,lwd=3)	## median
lines(xpred,ci[3,],col=3,lty=2)	## upper CI
lines(xpred,pi[1,],col=4,lty=2)	## lower PI
lines(xpred,pi[2,],col=4,lty=2)	## upper PI
