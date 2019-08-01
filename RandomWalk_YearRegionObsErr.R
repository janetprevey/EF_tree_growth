library(dplyr)
library(ecoforecastR)

precip <- read.csv("buckhorn_daily_precip.csv")
temps <- read.csv("buckhorn_daily_temps.csv")
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

RandomWalk_YearRegionObsErr <- "
model{
tau_obs ~ dgamma(a_obs,r_obs)
tau_r ~ dgamma(r1,r2)
tau_t ~ dgamma(t1,t2)

for(r in 1:NR){
alpha_r[r] ~ dnorm(0,tau_r)
x[1,r] ~ dnorm(x_ic,tau_ic)
}

for(t in 1:NT){
alpha_t[t] ~ dnorm(0,tau_t)
tau_add[t] ~ dgamma(a_add,r_add)
}

for(t in 2:NT){ ## t = year
for(r in 1:NR){
Ex[t,r] <- x[t-1,r] + alpha_r[r] + alpha_t[t]
x[t,r] ~ dnorm(Ex[t,r],tau_add[t])  ## process model
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

data <- list(region=region,NR=12,n=length(y),time=time,y=y,NT=10,x_ic=5,tau_ic=100,a_obs=1,r_obs=1,a_add=1,r_add=1,r1=0.1,r2=0.1,t1=0.1,t2=0.1)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_obs=5/var(y.samp))
}

j.model   <- jags.model(file = textConnection(RandomWalk_YearRegionObsErr),
                        data = data,
                        inits = init,
                        n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","alpha_r","alpha_t","tau_obs","tau_add"),
                            n.iter = 5000)

##evaluate model
#plot(jags.out)
#gelman.plot(jags.out)

##remove burnin
burnin = 500                                ## determine convergence
jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
#gelman.plot(jags.burn)

#acfplot(jags.burn)
effectiveSize(jags.burn)
summary(jags.burn)
out <- as.matrix(jags.burn)