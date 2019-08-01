library(dplyr)
library(ecoforecastR)

#read in tree growth data
data <- read.csv("buckhorn_growth_2007_2018.csv") %>%
  mutate(time = rep(1:10, 970))

#remove non-integer values from column "row"
sub <- data[data$row %%1 == 0,] 

#filter trees to see tag numbers of trees labeled "NP" (Not Planted)
NP <- filter(sub, alive == "NP") %>%
  select(tag_num)

# remove trees labled "NP" (Not Planted) in 2009 from all years, 
# remove trees that died during experiment, remove unneeded columns 
growth <- sub %>%
  arrange(plot) %>%
  filter(!tag_num %in% c(4909,56,4856,4887,2142,3682)) %>%
  filter(alive == "y") %>%
  select(-site, -site_num, -dead_year, -notes) 
  
# create separate mortality dataframe, where:
  # 2 = alive
  # 1 = dead
# This is to incorporate mortality into future models 
live_dead <- growth %>%
   mutate(mortality = as.integer(alive)) %>%
   select(plot, tag_num, year, mortality) 

RandomWalk_Regions <- "
model{
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_r ~ dgamma(r1,r2)
  tau_add ~ dgamma(a_add,r_add)
  
  for(r in 1:NR){
    alpha_r[r] ~ dnorm(0,tau_r)
    x[1,r] ~ dnorm(x_ic,tau_ic)
  }
  
  for(t in 2:NT){ ## t = year
    for(r in 1:NR){
      Ex[t,r] <- x[t-1,r] + alpha_r[r]
      x[t,r] ~ dnorm(Ex[t,r],tau_add)  ## process model
    }
  }
  
  for(i in 1:n){ ## i = individual trees
    y[i] ~ dnorm(x[time[i],region[i]],tau_obs)		        ## data model
  }
}
"
y <- growth$HT_cm
time <- growth$time
region <- as.integer(growth$region)

data <- list(region=region,NR=12,n=length(y),time=time,y=y,NT=10,x_ic=5,tau_ic=100,a_obs=1,r_obs=1,a_add=1,r_add=1,r1=0.1,r2=0.1)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)),tau_obs=5/var(y.samp))
}

j.model   <- jags.model(file = textConnection(RandomWalk_Regions),
                        data = data,
                        inits = init,
                        n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","alpha_r","tau_obs","tau_add"),
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
head(out)
nsamp <- 5000
samp <- sample.int(nrow(out),nsamp)
xpred <- seq(1,10,length=10)  					## sequence of x values we're going to
npred <- length(xpred)				##      make predictions for


alpha_matrix <- out[,1:12]
head(alpha_matrix)

intervals_by_region <- list()

for(i in seq_len(12)){# region loop
  ypred <- matrix(0.0,nrow=nsamp,ncol=npred)	## storage for predictive interval
  ycred <- matrix(0.0,nrow=nsamp,ncol=npred)	## storage for credible interval
  for(g in seq_len(nsamp)){ # sample loop
    theta <- out[samp[g],]
    alpha_i <- theta[i] 
    ycred[g,] <- theta[grep(paste(",",i,"]",sep = ""),names(theta))]+ alpha_i
    ypred[g,] <- rnorm(npred, ycred[g,], 1/sqrt(theta["tau_obs"]))
  }
  
  ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))  ## credible interval and median
  pi <- apply(ypred,2,quantile,c(0.025,0.975))		## prediction interval
  intervals_by_region[[i]]<- rbind(ci,pi)
}

#'
#' # Plot for each region over time
#'
## states  

head(intervals_by_region)
dev.off()
plot(seq(from=0,to=1200,length.out = 10)~xpred,type="n",ylab = "Height (cm)",xlab="Year")

for(s in 1:12){
  ciEnvelope(x = xpred,ylo = intervals_by_region[[s]][1,],
             yhi = intervals_by_region[[s]][3,],col=col.alpha(s,0.6))
  points(xpred,intervals_by_region[[s]][2,],col=s,pch=19)
  points(xpred,intervals_by_region[[s]][2,])
}

