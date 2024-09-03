
data {
  
  int n;
  vector[3] prizes;
  matrix[n,3] probLeft;
  matrix[n,3] probRight;
  int choiceLeft[n];
  
  int nalpha;
  vector[nalpha] alpha;
  
  vector[2] prior_r;
  vector[2] prior_lambda;
  
  int UseData;
  int UsePrior;
  
}

transformed data {
  
  matrix[n,3] dprob = probLeft-probRight;
  
}

parameters {
  
  // non-centered parameterization
  real theta_r;
  real theta_lambda;
  
}

transformed parameters {
  
  real r = prior_r[1]+prior_r[2]*theta_r;
  real lambda = exp(prior_lambda[1]+prior_lambda[2]*theta_lambda);
  
}

model {
  
  // prior
  if (UsePrior==1) {
    
    theta_r ~ std_normal();
    theta_lambda ~ std_normal();
  }
  
  // likelihood
  vector[3] u = pow(prizes,1.0-r)/(1.0-r);
  vector[3] unormalized = u/(u[3]-u[1]);
  
  if (UseData==1) {
    choiceLeft ~ bernoulli_logit(lambda*dprob*unormalized);
    
  }
  
  
  
}

generated quantities {
  
  vector[nalpha] CE = pow(alpha*pow(prizes[3],1.0-r)+(1.0-alpha)*pow(prizes[1],1.0-r),1.0/(1.0-r));
  
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
