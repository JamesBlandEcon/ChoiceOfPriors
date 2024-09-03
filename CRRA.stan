
data {
  
  int n; // number of observations
  vector[3] prizes; // the 3 monetary prizes ($10, $30, and $50 for HN2016). MUST BE ORDERED FROM SMALLEST TO LARGEST FOR THE CONTEXTUAL UTILITY NORMALIZATION TO WORK
  matrix[n,3] probLeft; // probabilities for the left lottery. Eah row is a different lottery
  matrix[n,3] probRight; // probabilities for the right lottery. Eah row is a different lottery
  int choiceLeft[n];// indicator variable for whether the left lottery was chosen
  
  /*
  Certainty equivalents:
  This code will evaluate an arbitrary number of certainty equivalents of the form:
      p x(highest prize),   (1-p)x(lowest prize)
  */
  int nalpha; // number of certainty equivalents to evaluate
  vector[nalpha] alpha; // the probability on the high prize to evaluate
  
  vector[2] prior_r; // prior for r, normal(mean,sd)
  vector[2] prior_lambda; // prior for lambda, lognormal(mean,sd)
  
  int UseData; // indicator variable. Set =1 to use the data, set =0 to ignore the likelihood and get draws from the prior
  int UsePrior; // indicator variable. Set =1 to use the prior, set =0 to ignore the prior (needed for maximum likelihood estimation)
  
}

transformed data {
  
  /* Expected utility maximizers in this setting only care about differences
  in probabilities. We can pre-compute these so we don't have to keep subtracting 
  one from the other for every iteration.
  */
  
  matrix[n,3] dprob = probLeft-probRight;
  
}

parameters {
  
  /* In an attempt to reduce the number of divergent transitions, I explored using 
  a non-centered parameterization. This is probably unnecessary, and didn't 
  fix the problem appreciably. 
  */
  
  // non-centered parameterization
  real theta_r;
  real theta_lambda;
  
}

transformed parameters {
  
  // converting the non-centered parameterization into the parameters we care about
  real r = prior_r[1]+prior_r[2]*theta_r;
  real lambda = exp(prior_lambda[1]+prior_lambda[2]*theta_lambda);
  
}

model {
  
  // prior
  if (UsePrior==1) {
    
    theta_r ~ std_normal();
    theta_lambda ~ std_normal();
  }
  
  // calculate the lilekihood
  vector[3] u = pow(prizes,1.0-r)/(1.0-r); // utility of each prize
  vector[3] unormalized = u/(u[3]-u[1]); // Contextual utility normalization
  
  if (UseData==1) {
    choiceLeft ~ bernoulli_logit(lambda*dprob*unormalized); // likelihood comtribution
    
  }
  
  
  
}

generated quantities {
  
  // Certainty equivalent
  vector[nalpha] CE = pow(alpha*pow(prizes[3],1.0-r)+(1.0-alpha)*pow(prizes[1],1.0-r),1.0/(1.0-r));
  
  /* This program also calculates the certainty equivalent of participating in 
  whole HN2016 experiment. That is, if someone were to make the 80 decisions 
  with behavior described parameters r and lambda, how much would they value the 
  experiment?
  
  */
  real CEexperiment;
  
  { 
    vector[3] u = pow(prizes,1.0-r)/(1.0-r);
    vector[3] unormalized = u/(u[3]-u[1]);
    vector[n] pL = inv_logit(lambda*dprob*unormalized);
    
    row_vector[3] ExperimentLottery;
    for (xx in 1:3) {
      ExperimentLottery[xx] = mean(pL.*probLeft[,xx]+(1.0-pL).*probRight[,xx]);
    }
    
    CEexperiment = pow(ExperimentLottery*pow(prizes,1.0-r),1.0/(1.0-r));
    
    
  }
  
  
}
