
library(tidyverse)
library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

source("shortcode/loadData.R")
  
model<-"shortcode/CRRA.stan" |> 
  stan_model()
  
  
prior_r<-c(0.27,0.36)
prior_lambda<-c(3.45,0.59)

SimSize<-1000

ITER<-10000
ADAPT_DELTA <- 0.8
CORES <- 4

alpha<-c(0.2,0.5,0.8)

reRunSimulation<-TRUE
ReRunEstimates<-FALSE


################################################################################
# Draws from the priors
################################################################################

file<-"shortoutputs/priors.rds"

if (!file.exists(file)) {
  
  PRIORS<-tibble()
  
  d <- D |>
    filter(id==ii)
  
  
  # Estimate model with the calibrated prior
  dStan<-list(
    n = dim(d)[1],
    prizes = c(10,30,50),
    probLeft = cbind(d$probL1,d$probL2,d$probL3),
    probRight = cbind(d$probR1,d$probR2,d$probR3),
    choiceLeft = d$choiceLeft,
    
    nalpha = length(alpha),
    alpha = alpha,
    
    prior_r = prior_r,
    prior_lambda = prior_lambda,
    
    UseData=0,
    UsePrior=1
  )
  
  Fit<-model |>
    sampling(
      data=dStan,seed=42,
      cores=1
    )
  
  addThis<-tibble(
    r = extract(Fit)$r,
    lambda = extract(Fit)$lambda,
    CE1 = extract(Fit)$CE[,1],
    CE2 = extract(Fit)$CE[,2],
    CE3 = extract(Fit)$CE[,3],
    CEexperiment = extract(Fit)$CEexperiment,
    prior = "calibrated"
  )
  
  PRIORS<-rbind(PRIORS,addThis)
  
  dStan$prior_r[2]<-dStan$prior_r[2]*priorFactor
  dStan$prior_lambda[2]<-dStan$prior_lambda[2]*priorFactor
  
  Fit<-model |>
    sampling(
      data=dStan,seed=42,
      cores=1
    )
  addThis<-tibble(
    r = extract(Fit)$r,
    lambda = extract(Fit)$lambda,
    CE1 = extract(Fit)$CE[,1],
    CE2 = extract(Fit)$CE[,2],
    CE3 = extract(Fit)$CE[,3],
    CEexperiment = extract(Fit)$CEexperiment,
    prior = "inflated"
  )
  
  PRIORS<-rbind(PRIORS,addThis)
  
  PRIORS |> saveRDS(file)
} else {
  print("Skipping prior simulation: file already exists")
  
  PRIORS<-file |>
    readRDS()
}


set.seed(42)

###############################################################################
# Estimate EUT model from data in HN2016.
# I do this for the calibrated prior and the inflated prior
###############################################################################

file<-"shortoutputs/estimates.rds"


if (!file.exists(file) | ReRunEstimates) {
  ESTIMATES<-tibble()
  
  for (ii in unique(D$id)) {
    
  print(paste("Estimating for id =",ii,"--------------------------------"))
    
    
    d <- D |>
      filter(id==ii)
    
    
    # Estimate model with the calibrated prior
    dStan<-list(
      n = dim(d)[1],
      prizes = c(10,30,50),
      probLeft = cbind(d$probL1,d$probL2,d$probL3),
      probRight = cbind(d$probR1,d$probR2,d$probR3),
      choiceLeft = d$choiceLeft,
      
      nalpha = length(alpha),
      alpha = alpha,
      
      prior_r = prior_r,
      prior_lambda = prior_lambda,
      
      UseData=1,
      UsePrior=1
    )
      
    Fit<-model |>
        sampling(
          data=dStan,seed=42,
          # This runs faster without the overhead from parallelization
          cores=1,
          control = list(adapt_delta = 0.99)
        )
    
    divergences<-sum(Fit |> get_divergent_iterations())
    
    addThis<-summary(Fit)$summary |>
      data.frame() |>
      rownames_to_column(var="par") |>
      mutate(
        id = ii,
        prior = "calibrated",
        divergences = divergences
      )
    
    ESTIMATES<-rbind(
      ESTIMATES,addThis
    )
    
    # Estimate the model with the inflated prior for both parameters
    
    dStan$prior_r[2]<-prior_r[2]*priorFactor
    dStan$prior_lambda[2]<-prior_lambda[2]*priorFactor
    
    Fit<-model |>
      sampling(
        data=dStan,seed=42,
        # This runs faster without the overhead from parallelization
        cores=1,
        control = list(adapt_delta = 0.99)
      )
    
    divergences<-sum(Fit |> get_divergent_iterations())
    
    addThis<-summary(Fit)$summary |>
      data.frame() |>
      rownames_to_column(var="par") |>
      mutate(
        id = ii,
        prior = "both inflated",
        divergences = divergences
      )
    
    ESTIMATES<-rbind(
      ESTIMATES,addThis
    )
    
    # Estimate the model with the inflated prior for r
    
    dStan$prior_r[2]<-prior_r[2]*priorFactor
    dStan$prior_lambda[2]<-prior_lambda[2]
    
    Fit<-model |>
      sampling(
        data=dStan,seed=42,
        # This runs faster without the overhead from parallelization
        cores=1,
        control = list(adapt_delta = 0.99)
      )
    
    divergences<-sum(Fit |> get_divergent_iterations())
    
    addThis<-summary(Fit)$summary |>
      data.frame() |>
      rownames_to_column(var="par") |>
      mutate(
        id = ii,
        prior = "r inflated",
        divergences = divergences
      )
    
    ESTIMATES<-rbind(
      ESTIMATES,addThis
    )
    
    # Estimate the model with the inflated prior for LAMBDA
    
    dStan$prior_r[2]<-prior_r[2]
    dStan$prior_lambda[2]<-prior_lambda[2]*priorFactor
    
    Fit<-model |>
      sampling(
        data=dStan,seed=42,
        # This runs faster without the overhead from parallelization
        cores=1,
        control = list(adapt_delta = 0.99)
      )
    
    divergences<-sum(Fit |> get_divergent_iterations())
    
    addThis<-summary(Fit)$summary |>
      data.frame() |>
      rownames_to_column(var="par") |>
      mutate(
        id = ii,
        prior = "lambda inflated",
        divergences = divergences
      )
    
    ESTIMATES<-rbind(
      ESTIMATES,addThis
    )
      
  }
  
  
  ESTIMATES |>
    saveRDS(file)
  
  
} else {
  print("Skipping individual-level estimation: file already exists")
  
}


###############################################################################
# Simulate data from the calibrated prior and repeat this exercise
###############################################################################

file<-"shortoutputs/simulation.rds"


if (!file.exists(file)| reRunSimulation) {
  
  SIMULATION<-tibble()

  d<-D |> filter(id==1)
  
  dprob<-cbind(d$probL1,d$probL2,d$probL3)-cbind(d$probR1,d$probR2,d$probR3)
  prizes<-c(10,30,50)
  
  inv.logit<-function(x) 1/(1+exp(-x))
  
  simData<-function(r,lambda) {
    
    u<-(prizes^(1-r))/(1-r)
    u<-u/(u[3]-u[1])
    
    pL<-inv.logit(
      lambda*dprob%*%u
    )
    
    (1*(runif(dim(dprob)[1])<pL)) |> as.vector()
  }
  
  priors<-PRIORS |> filter(prior=="calibrated") |> select(-prior)
  
  
  good_runs<-0
  
  ss<-0
  
  while (good_runs < SimSize) {
    
    ss<-ss+1
    print(paste("Simulation step",ss,"of",SimSize))
  
    
    
    truth<-c(NA,NA,priors[ss,] |> as.vector() |> unlist() |> c(NA))
    
    r<-truth["r"]
    lambda<-truth["lambda"]
    
    # Estimate model with the calibrated prior
    dStan<-list(
      n = dim(d)[1],
      prizes = c(10,30,50),
      probLeft = cbind(d$probL1,d$probL2,d$probL3),
      probRight = cbind(d$probR1,d$probR2,d$probR3),
      choiceLeft = simData(r,lambda),
      
      nalpha = length(alpha),
      alpha = alpha,
      
      prior_r = prior_r,
      prior_lambda = prior_lambda,
      
      UseData=1,
      UsePrior=1
    )
    
    Fit<-model |>
      sampling(
        data=dStan,seed=42,
        cores=CORES,
        control = list(adapt_delta = ADAPT_DELTA),
        iter = ITER,
      )
    
    divergences<-sum(Fit |> get_divergent_iterations())
    
    addThis<-summary(Fit)$summary |>
      data.frame() |>
      rownames_to_column(var="par") |>
      mutate(
        id = ss,
        prior = "calibrated",
        divergences = divergences,
        truth=truth
      )
    
    SIMULATION<-rbind(
      SIMULATION,addThis
    )
    
    # Estimate the model with the 2x inflated prior
    
    dStan$prior_r[2]<-prior_r[2]*2
    dStan$prior_lambda[2]<-prior_lambda[2]*2
    
    Fit<-model |>
      sampling(
        data=dStan,seed=42,
        cores=CORES,
        control = list(adapt_delta = ADAPT_DELTA),
        iter = ITER,
      )
    
    divergences<-sum(Fit |> get_divergent_iterations())
    
    addThis<-summary(Fit)$summary |>
      data.frame() |>
      rownames_to_column(var="par") |>
      mutate(
        id = ss,
        prior = "x2",
        divergences = divergences,
        truth=truth
      )
    
    SIMULATION<-rbind(
      SIMULATION,addThis
    )
    
    # Estimate the model with the 8x inflated prior
    
    dStan$prior_r[2]<-prior_r[2]*8
    dStan$prior_lambda[2]<-prior_lambda[2]*8
    
    Fit<-model |>
      sampling(
        data=dStan,seed=42,
        cores=CORES,
        control = list(adapt_delta = ADAPT_DELTA),
        iter = ITER,
      )
    
    divergences<-sum(Fit |> get_divergent_iterations())
    
    addThis<-summary(Fit)$summary |>
      data.frame() |>
      rownames_to_column(var="par") |>
      mutate(
        id = ss,
        prior = "x8",
        divergences = divergences,
        truth=truth
      )
    
    SIMULATION<-rbind(
      SIMULATION,addThis
    )
    
    
    # MLE -------------------------------------------------------------------
    
    dStan$UsePrior<-0
    MLE<-model |>
      optimizing(data=dStan,hessian = TRUE,
                 init = list(r = truth["r"], lambda = truth["lambda"])
                 )
    
    addThis<-tibble(
      par = names(MLE$par),
      mean = MLE$par,
      se_mean = NA,
      sd = NA, #c(NA,NA,-MLE$hessian |> solve() |> sqrt() |> diag(), rep(NA,4)),
      X2.5. = NA,
      X25. = NA,
      X50. = NA,
      X75.  = NA,
      X97.5. = NA,
      n_eff = 10000,
      Rhat = 1,
      id = ss,
      prior = "MLE",
      divergences = MLE$return_code,
      truth = truth[1:8]
    )
    
    SIMULATION<-rbind(SIMULATION,addThis)
    
    
    check <- SIMULATION |>
      group_by(id) |>
      mutate(
        drop_divergent = max(divergences)!=0,
        drop_Rhat = max(Rhat)>1.01,
        drop_n_eff = min(n_eff)<1000,
        drop = max(drop_divergent+drop_Rhat+drop_n_eff)>=1
      ) |>
      ungroup() |>
      filter(prior=="calibrated",par=="r") |>
      summarize(dropped = sum(drop),
                divergent = sum(drop_divergent),
                Rhat = sum(drop_Rhat),
                n_eff = sum(drop_n_eff),
                n = n()
      )
    
    check |> print()
    
    good_runs<-check$n-check$dropped
    
    
    SIMULATION |>
      saveRDS(file)
    
  }
  
  
  

} else {
  print("skipping simulation: file already exists")
}

