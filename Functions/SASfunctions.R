# All the SAS functions

# Power law function
fSAS_pl = function(Ps, par, ...) {
  # Compute power law function
  # Ps = cumalative storage age distribution
  # pars = list with parameter k (exponent)
  
  k = par$k
  
  # Constraints
  if(k<0) stop('pars must be positive')
  
  # Compute omega
  Om = Ps^k
  
  return(Om)
}

# Compute power law function depending on system state wi.
fSAS_pltv = function(Ps, par, wi) {
  # Ps = cumalative storage age distribution
  # pars = list with kmin and kmax
  # wi   = system state (e.g. deltaS scaled to 0 to 1)
  
  kmin=par$kmin
  kmax=par$kmax
  
  # Constraints
  if(any(c(kmin,kmax)<0)) stop('pars must be positive')
  
  k = kmin+(1-wi)*(kmax-kmin)
  
  # Compute omega
  Om = Ps^k
  
  return(Om)
}

# Cumalative beta distribution function
fSAS_beta = function(Ps, par, ...) {
  # Returns 
  # Ps = cumalative storage age distribution
  # par = list of parameters alpha and beta
  
  # Get shape pars
  a = par$alpha
  b = par$beta
  
  # Compute omega
  Om = pbeta(Ps, a, b)
  return(Om)
}
