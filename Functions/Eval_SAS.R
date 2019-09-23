# Evaluate SAS performance

evalSAS = function(obs, sim) {
  # Obs: observed
  # Sim: simulated
  # Returns: Nash-Sutcliffe Efficiency
  #          Kling-Gupta Efficiency
  #          Lin's concordance
  
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
  
  # Out
  out = list(NSE  = NSE,
             KGE  = KGE,
             LCCC = LCCC)
  return(out)
}