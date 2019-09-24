# Evaluate SAS performance

evalSAS = function(obs, sim, npar=0) {
  # obs: observed
  # sim: simulated
  # npar: supply amount of pars optimised for AIC calculation
  # Returns: Nash-Sutcliffe Efficiency
  #          Kling-Gupta Efficiency
  #          Lin's concordance
  #          AIC
  
  # Check for NA's in obs data
  if(sum(is.na(obs))>0) {
    sim = sim[!is.na(obs)]
    obs = obs[!is.na(obs)]
  }
  
  # NSE
  NSE = 1- mean((obs - sim)^2) / var(obs)
  
  # KGE
  KGE <- 1 - sqrt( (cor(sim, obs) - 1)^2 +
                   (mean(sim)/mean(obs) - 1)^2 +
                   (sd(sim)/sd(obs) - 1)^2 )
  # LCCC
  n     = length(obs)
  sx2   = var(sim) * (n-1) / n      # Variance sim
  sy2   = var(obs) * (n-1) / n      # Variance obs
  sxy   = cov(sim, obs) * (n-1) / n # Covariance
  LCCC  = 2 * sxy / (sx2 + sy2 + (mean(sim)-mean(obs))^2)
  
  # RMSE
  RMSE = sqrt( mean( (obs-sim)^2 ) )
  
  # AIC
  if(npar>0) {
    AIC  = length(obs)*log(RMSE) + 2*npar

    out = list(NSE  = NSE,
               KGE  = KGE,
               LCCC = LCCC,
               RMSE = RMSE,
               AIC  = AIC)
  } else {
    # Out
    out = list(NSE  = NSE,
               KGE  = KGE,
               LCCC = LCCC,
               RMSE = RMSE)
  }
  
  return(out)
}


