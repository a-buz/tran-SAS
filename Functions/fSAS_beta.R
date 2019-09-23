fSAS_beta = function(Ps, par) {
  # Returns cumalative beta distribution function
  # Ps = cumalative storage age distribution
  # par = list of parameters alpha and beta
  
  # Get shape pars
  a = par$alpha
  b = par$beta
  
  # Compute omega
  Om = pbeta(Ps, a, b)
  return(Om)
}
