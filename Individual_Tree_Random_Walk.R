#Load libraries and data
library(dplyr)
library(rjags)
library(ecoforecastR)

datadir <- "/Users/rzabramoff/Documents/Forecasting/EF_tree_growth/"
precip <- read.csv(paste0(datadir,"buckhorn_daily_precip.csv"))
temps <- read.csv(paste0(datadir,"buckhorn_daily_temps.csv"))
growth <- read.csv(paste0(datadir,"buckhorn_growth_2007_2018.csv"))


#Fit a random walk to one tree
RandomWalk = "
model{

#### Data Model
for(t in 1:n){
y[t] ~ dnorm(x[t],tau_obs)
}

#### Process Model
for(t in 2:n){
x[t]~dnorm(x[t-1],tau_add)
}

#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
}
"

y <- growth[growth$tag_num==as.character(unique(growth$tag_num)[1]),]$HT_cm
time <- growth[growth$tag_num==as.character(unique(growth$tag_num)[1]),]$year

data <- list(y=y,n=length(y),x_ic=5,tau_ic=100,a_obs=1,r_obs=1,a_add=1,r_add=1)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)),tau_obs=5/var(y.samp))
}

j.model   <- jags.model(file = textConnection(RandomWalk),
                        data = data,
                        inits = init,
                        n.chains = 3)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
plot(jags.out)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)

burnin = 2000                               ## determine convergence
jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
plot(jags.burn)                             ## check diagnostics post burn-in

time.rng = c(1,length(time)) ## adjust to zoom in and out
out <- as.matrix(jags.burn)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y, na.rm=T),ylab="Height (cm)",xlim=time[time.rng])
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)
points(time, ci[2,], col=2)

plot(time,ci[2,],type='n',ylim=range(ci),ylab="Height (cm)",xlim=time[time.rng])
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)
points(time, ci[2,], col=2)
