# The circglmbayes package needs to be installed from github
# install.packages("devtools")
# devtools::install_github("keesmulder/circglmbayes")
# 
# For brms on windows we need to install the latest development version of rstan,
# i.e. not the one from CRAN. This is not needed on MacOS.
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(circglmbayes)
library(coda) #for mcmc tools to use with circglmbayes
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)

par(mar = c(1, 1, 1, 1))

setwd("C:/Users/skly5321/Downloads/ChapterOne/ChapterOne/Data")


# The link function for a circular GLM is the tan half function:
#   tan_half(x) = tan(x/2)
# So the inverse link function is the inverse tan half

inv_tan_half <- function(x) {
  return(2 * atan(x))
}

# Conversion functions for day of year to circular and vice versa. These convert
# day of year to the range -pi:pi, where -pi and pi are the same point on the
# circle. This does not account for leap years (for leap years, day of year
# needs to be calculated to include 29 Feb and divided by 366 instead of 365).

day2circ <- function(day_of_year) {
  return(2 * pi * day_of_year / 365 - pi)
}

circ2day <- function(circ) {
  return(365 * (circ + pi) / (2 * pi))
}

# There are at least two versions of circular GLMs that use the tan half link
# function in different ways. The most sensible seems to be the version
# described by Mulder and Klugkist 2017 and implemented in the `circglmbayes`
# package. While described as a "GLM", this is not really a GLM (so it can't be
# fit in any existing GLM packages that allow for custom GLMs) and is instead a
# nonlinear model. In this model, the intercept terms are **outside** the link
# function. The key reason this model works the best is that it wraps properly
# on the circular scale. For a circular response y on -pi:pi and a single
# continuous predictor x, the model is:

# $$
# y = \beta_0 + inv_tan_half(\beta_1 * x)
# $$ 

# Let's see how this model works.


v <- c("Amazon", "Bangladesh", "Bolivia", "Caatinga", "Cameroon", "CochaCashu(new)", "FrenchPolynesia", "Ghana", "Guinnea", "JatunSacha", "LasCruces", "MadhyaPradesh", "Miombo", "Myanmar", "Udzungwa", "WesternGhatsIndia")

traceplots_list <- list()
finalplot_list <- list()
coef_list <- list()

for (i in 1:16){
  df <- read.csv(paste(v[i],".csv", sep = ""))
  df$species <- factor(df$species)
  spec <- levels(df$species)
  for(j in length(spec)){
    plant2 <- na.omit(df[df$species == spec[j],])
    year <- plant2$Year
    year_s <- scale(year) #year scaled and centered
    year_center <- attr(year_s, "scaled:center")
    year_scale <- attr(year_s, "scaled:scale")
    circ <- day2circ(plant2$Day.of.the.year)
    sdat <- data.frame(year, year_s, circ)
    nchains <- 4
    chains <- list()
    for (k in 1:nchains ) {
      fit <- circGLM(circ ~ year_s, data=sdat)
      chains[[k]] <- fit$all_chains
    }
    chains <- mcmc.list(chains)
    traceplots_list[[spec[j]]] <- traceplot(chains)
    autocorr.plot(chains, ask=FALSE)
    gelman.diag(chains[,"b0_chain"])     
    gelman.diag(chains[,"kp_chain"])     
    gelman.diag(chains[,"bt_chain"])  
    fit <- circGLM(circ ~ year_s, data=sdat, burnin=200, thin=30, Q=2500)
    fit
    samples <- fit$all_chains
    samplesdf <- data.frame(samples)
    names(samplesdf) <- names(samplesdf) %>% 
      gsub("_chain", "", .)
    samplesdf %>% 
      select("b0","kp","bt") %>% 
      pivot_longer(cols=everything(), names_to="parameter", values_to="sample_value") %>% 
      ggplot() +
      geom_histogram(mapping=aes(x=sample_value, y=after_stat(density)), bins=75) +
      facet_wrap(vars(parameter), scales="free")
    
    mean(samplesdf$b0)
    mean(samplesdf$kp)
    mean(samplesdf$bt)
    
    hpdi <- function (samp, prob = 0.95) {
      vals <- sort(samp)
      nsamp <- length(vals)
      gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
      init <- 1:(nsamp - gap)
      inds <- which.min(vals[init + gap,drop=FALSE] - vals[init, drop=FALSE])
      ans <- cbind(lower=vals[inds], upper=vals[inds + gap])
      return(ans)
    }
    
    hpdi(samplesdf$b0, prob=0.95)
    hpdi(samplesdf$kp, prob=0.95)
    hpdi(samplesdf$bt, prob=0.95)
    
    quantile(samplesdf$b0, prob=c(0.025,0.975))
    quantile(samplesdf$kp, prob=c(0.025,0.975))
    quantile(samplesdf$bt, prob=c(0.025,0.975))
    
    
    year <- seq(from=min(sdat$year), to=max(sdat$year), by=1) 
    year_s <- scale(year, center=year_center, scale=year_scale)
    n <- length(year)
    results <- matrix(NA, nrow=n, ncol=5) 
    colnames(results) <- c("mnmu","mulo95","muhi95","ppdlo95","ppdhi95")
    
    for ( i in 1:n ) {
      
      mu <- samplesdf$b0 + inv_tan_half(samplesdf$bt * year_s[i])
    
      ppd <- rvon_mises(n=length(mu), mu=mu, kappa=samplesdf$kp)
      
      results[i,1] <- mean(mu)
      #results[i,2:3] <- hpdi(mu, prob=0.95)
      results[i,2:3] <- quantile(mu, prob=c(0.025,0.975)) #CPI
      #results[i,4:5] <- hpdi(ppd, prob=0.95)
      results[i,4:5] <- quantile(ppd, prob=c(0.025,0.975)) #CPI
    }
    results <- circ2day(results) #transform to day of year
    preds <- data.frame(year, year_s, results)
    rm(year, year_s, n, results, mu, ppd) #clean up
    
    coef_list[[spec[j]]] <- coef(fit)
    
    finalplot_list[[spec[j]]] <- preds %>%
      ggplot() +
      geom_ribbon(mapping=aes(x=year, ymin=mulo95, ymax=muhi95), alpha=0.2) +
      geom_point(data=sdat, 
                 mapping=aes(x=year, y=plant2$Day.of.the.year)) +
      geom_line(mapping=aes(x=year, y=mnmu)) +
      geom_line(mapping=aes(x=year, y=ppdlo95), lty=2) +
      geom_line(mapping=aes(x=year, y=ppdhi95), lty=2)
  }
}

















setwd("C:/Users/skly5321/Downloads/ChapterOne/ChapterOne/Data")
plant2 <- read.csv("leptochloa_panicoides.csv")

plot(plant2$Day.of.the.year,plant2$Year)




# If we were to examine a histogram of the circular-transformed data (i.e.
# -pi:pi scale), it would look like this:


# Package in a data frame

year <- plant2$Year
year_s <- scale(year) #year scaled and centered
year_center <- attr(year_s, "scaled:center")
year_scale <- attr(year_s, "scaled:scale")
circ <- day2circ(plant2$Day.of.the.year)


sdat <- data.frame(year, year_s, circ)

plot(year, circ)

# Fit the model with `circGLM()` from the `circglmbayes` package. The training
# algorithm is a mix of Gibbs sampling, fast rejection sampling, and
# Metropolis-Hastings. By default, this function uses 1 chain and has no burn in
# period. There are also no convergence diagnostics by default. So we will build
# these analysis steps in. We'll run multiple chains, keep the results in a
# list, and use the diagnostic tools from the `coda` package since this function
# was designed to work with those.

nchains <- 4
chains <- list()
for (i in 1:nchains ) {
  fit <- circGLM(circ ~ year_s, data=sdat)
  chains[[i]] <- fit$all_chains
}
chains <- mcmc.list(chains) #convert to a mcmc.list, a `coda` data structure

# The traceplots help us to see when the burn in period is over. It looks like a
# burnin period of 200 would be ample.

traceplot(chains)

# Autocorrelation tells us how much to thin the samples so that we have samples
# that are effectively independent. It looks like the autocorrelation has
# completely gone away after about 20-30 samples.

par(mar = c(1, 1, 1, 1))

autocorr.plot(chains, ask=FALSE)

# Gelman's R-hat diagnostic. These are all close to 1, indicating convergence.

gelman.diag(chains[,"b0_chain"])     
gelman.diag(chains[,"kp_chain"])     
gelman.diag(chains[,"bt_chain"])     

# Based on the diagnostics, we'll allow a burnin of 200, thin by 30 and collect
# 2500 samples, all in a single chain.

fit <- circGLM(circ ~ year_s, data=sdat, burnin=200, thin=30, Q=2500)

# The fitted model recovers the known parameters quite well

fit

# Extract the samples and package in a data frame with nicer parameter names
samples <- fit$all_chains
samplesdf <- data.frame(samples)
names(samplesdf) <- names(samplesdf) %>% 
  gsub("_chain", "", .)

# Use the samples to form inferences. This code is largely similar to the code
# from 10_8_ants_bayesian_GLM.

# Plot the posterior distributions for the parameters. bt is our original b1
# while kp is kappa.

samplesdf %>% 
  select("b0","kp","bt") %>% 
  pivot_longer(cols=everything(), names_to="parameter", values_to="sample_value") %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=sample_value, y=after_stat(density)), bins=75) +
  facet_wrap(vars(parameter), scales="free")

# Posterior means

mean(samplesdf$b0)
mean(samplesdf$kp)
mean(samplesdf$bt)

# Posterior HPDIs

hpdi <- function (samp, prob = 0.95) {
  vals <- sort(samp)
  nsamp <- length(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- which.min(vals[init + gap,drop=FALSE] - vals[init, drop=FALSE])
  ans <- cbind(lower=vals[inds], upper=vals[inds + gap])
  return(ans)
}

hpdi(samplesdf$b0, prob=0.95)
hpdi(samplesdf$kp, prob=0.95)
hpdi(samplesdf$bt, prob=0.95)

# Posterior CPIs (probably more stable than hpdi)

quantile(samplesdf$b0, prob=c(0.025,0.975))
quantile(samplesdf$kp, prob=c(0.025,0.975))
quantile(samplesdf$bt, prob=c(0.025,0.975))

# Regression and prediction intervals

# Initialize variables and storage
year <- seq(from=min(sdat$year), to=max(sdat$year), by=1) #range for year
year_s <- scale(year, center=year_center, scale=year_scale)
n <- length(year)
results <- matrix(NA, nrow=n, ncol=5) #to store hpdi values and mean
colnames(results) <- c("mnmu","mulo95","muhi95","ppdlo95","ppdhi95")

coef(fit)

# For each year, form the posterior
for ( i in 1:n ) {
  
  # Use inverse link for samples of the posterior \mu
  mu <- samplesdf$b0 + inv_tan_half(samplesdf$bt * year_s[i])
  
  # Sample from von mises to get the posterior predictive distribution
  ppd <- rvon_mises(n=length(mu), mu=mu, kappa=samplesdf$kp)
  
  # Mean and intervals of these samples
  results[i,1] <- mean(mu)
  #results[i,2:3] <- hpdi(mu, prob=0.95)
  results[i,2:3] <- quantile(mu, prob=c(0.025,0.975)) #CPI
  #results[i,4:5] <- hpdi(ppd, prob=0.95)
  results[i,4:5] <- quantile(ppd, prob=c(0.025,0.975)) #CPI
}
results <- circ2day(results) #transform to day of year
preds <- data.frame(year, year_s, results)
rm(year, year_s, n, results, mu, ppd) #clean up

preds %>%
  ggplot() +
  geom_ribbon(mapping=aes(x=year, ymin=mulo95, ymax=muhi95), alpha=0.2) +
  geom_point(data=sdat, 
             mapping=aes(x=year, y=plant2$Day.of.the.year)) +
  geom_line(mapping=aes(x=year, y=mnmu)) +
  geom_line(mapping=aes(x=year, y=ppdlo95), lty=2) +
  geom_line(mapping=aes(x=year, y=ppdhi95), lty=2)

# We can fit the same model with stan using the `brms` package. However it takes
# a lot longer, is somewhat fragile and needs a quite informative prior for
# \beta_0. We can probably make this work better with better priors. This is a
# nonlinear model, so needs to be set up as such in `brm()`. The parameter
# estimates are very similar.

prior_circ <- c(prior(normal(0, 0.02), nlpar = "b0"),
                prior(normal(0, 1), nlpar = "b1"))
fit_brms <- brm(bf(circ ~ b0 + 2 * atan(b1 * year_s), b0 ~ 1, b1 ~ 1, nl=TRUE),
                data=sdat, family=von_mises(link="identity"), prior=prior_circ, 
                iter=2000)
summary(fit_brms)