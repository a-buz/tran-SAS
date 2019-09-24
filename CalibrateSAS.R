# Simple example of SAS parameter calibration

# Libraries
library(tidyverse)

# Functions
source('Models/SAS_EFs.R')
source('Functions/SASfunctions.R')
source('Functions/Eval_SAS.R')

#------------------------------------------------------------------------------#
# Data
#------------------------------------------------------------------------------#
# Read in hourly test data that was generated from the original tran-SAS package
data = read_csv('Example.csv')
# Data has columns: Date, P, Q, ET, Cin, wi, Cout

# Set missing data to NA
data$Cout = ifelse(data$Cout==-999, NA, data$Cout)

# Aggregate hourly data to daily timestep
data = data %>% 
  mutate(Date = as.Date(Date, tz='UTC')) %>% 
  group_by(Date) %>% 
  summarise(Cin = weighted.mean(Cin, P), # Conserve mass input of Cin
            P  = sum(P),
            Q  = sum(Q),
            ET = sum(ET),
            wi = mean(wi),
            Cout = mean(Cout, na.rm=TRUE)) %>% 
  select(Date, P, Q, ET, Cin, wi, Cout)

# Fix tracer input for prettier plotting later on
getCin = data$Cin[!is.na(data$Cin)][1] # First non NA value
for(i in 1:nrow(data)) {
  if(data$P[i] > 0)  {
    getCin = data$Cin[i]
  } else {
    data$Cin[i] = getCin
  }
}

#------------------------------------------------------------------------------#
# Parameters that won't be calibrated
#------------------------------------------------------------------------------#
# Initial conditions
C_S0 = 0          # Initial concent of tracer in storage

# Other
f_thresh = 1 # Fraction of rank storage after which the storage is sampled uniformly

# Spinup
spinup = list(spinStart = '2012-09-30',   # Start of spinup
              spinEnd   = '2013-09-29',   # End of spinup
              spin_n    = 3)              # How many times to repeat, 0 equals no repeats

# Dates to sample age distributions
datesel = c('2015-08-15')

#------------------------------------------------------------------------------#
# Set up calibration function
#------------------------------------------------------------------------------#
calSAS = function(calPars, objFun='NSE') {
  # Parameters
  Qpars  = list(fun = 'fSAS_beta', pars = list(alpha = 1.0, beta = calPars[1]))
  ETpars = list(fun = 'fSAS_beta', pars = list(alpha = calPars[2], beta = 1.0))
  
  # Initial conditions
  S0   = calPars[3] # Storage amount
  
  # Combine parameters into a list
  pars = list(Q=Qpars, ET=ETpars, S0=S0, C_S0=C_S0, f_thresh=f_thresh,
              spinup=spinup, datesel=datesel)
  
  # Run model
  runMod = SAS_EFs(pars, data)
  
  # Return ObjFun
  #------------------------------------------------------------------------------#
  results   = evalSAS(obs = data$Cout, sim = runMod$C_Q, npar=length(calPars))
  objFunVal = results[objFun][[1]]
  
  return(objFunVal)
}

# Optimise the NSE
calResults = optim(par=c(0.3, 0.3, 750), # Initial values, Qbeta, ETalpha, S0
                   fn=calSAS,
                   method='L-BFGS-B',
                   lower=c(0,0,500), upper=c(1,1,100),
                   control=list(fnscale=-1, trace=1),
                   objFun='NSE')

#------------------------------------------------------------------------------#
# Results
#------------------------------------------------------------------------------#

# Get the results of the best calibration and rerun the model
Qpars  = list(alpha = 1.0, beta = calResults$par[1]) 
ETpars = list(alpha = calResults$par[2], beta = 1.0) 
S0   = calResults$par[3] 

pars = list(Q=Qpars, ET=ETpars, S0=S0, C_S0=C_S0, f_thresh=f_thresh,
            spinup=spinup, datesel=datesel)

runMod = SAS_EFs(pars, data)

# Matches ObjFun score
evalSAS(obs = data$Cout, sim = runMod$C_Q, npar=3) # npars optimised

# Plot
data$C_Q = runMod$C_Q
measC_Q  = ifelse(is.na(data$Cout), NA, data$C_Q)

ggplot(data, aes(x=Date)) +
  geom_line(aes(y=Cin, col='C in')) +
  geom_point(aes(y=Cout, col='C out')) +
  geom_line(aes(y=runMod$C_Q, col='Sim C_Q')) +
  geom_point(aes(y=measC_Q, col = 'Sim C_Q during meas')) +
  labs(col='', y='C') + 
  theme(legend.position = 'bottom')

# Plot obs vs sim
ggplot(data) + 
  geom_point(aes(Cout, C_Q)) +
  geom_abline(aes(slope=1, intercept=0)) +
  lims(x=c(45, 55), y=c(45, 55)) +
  labs(x='Obs', y='Sim')
