#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# R rewrite of tran-SAS package
# Alexander Buzacott

# Original MatLab version and citation:
# Benettin, P., & Bertuzzo, E. (2018). tran-SAS v1.0: a numerical model 
# to compute catchment-scale hydrologic transport using StorAge Selection 
# functions. Geoscientific Model Development, 11(4), 1627:1639. 
# https://doi.org/10.5194/gmd-11-1627-2018
# The codes implements the transport model described by:
# Rinaldo, A., Benettin, P., Harman, C. J., Hrachowitz, M., McGuire, K. J.,
# van der Velde, Y., Bertuzzo, E., and Botter, G. (2015). Storage selection 
# functions: A coherent framework for quantifying how catchments store and 
# release water and solutes. Water Resources Research, 51(6), 4840�4847. 
# http://doi.org/10.1002/2015WR017273
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Libraries
library(tidyverse) # To simplify aggregations and plotting

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
# Parameters
#------------------------------------------------------------------------------#
# Select functions for Q and ET and assign parameters
# Power law
# Qpars  = list(fun = 'fSAS_pl', pars = list(k=0.7))
# ETpars = list(fun = 'fSAS_pl', pars = list(k=0.7))

# Power law with system wetness
# Qpars  = list(fun = 'fSAS_pltv', pars = list(kmin=0.3, kmax=0.9))
# ETpars = list(fun = 'fSAS_pl',   pars = list(k=1))

# Beta distribution
Qpars  = list(fun = 'fSAS_beta', pars = list(alpha=1.0, beta=0.3))
ETpars = list(fun = 'fSAS_beta', pars = list(alpha=0.3, beta=1.0))

# Initial conditions
S0   = 1000 # Initial storage amount
C_S0 = 0    # Initial concent of tracer in storage

# Other
f_thresh = 1 # Fraction of rank storage after which the storage is sampled uniformly

# Spinup
spinup = list(spinStart = '2012-09-30',   # Start of spinup
              spinEnd   = '2013-09-29',   # End of spinup
              spin_n    = 3)              # How many times to repeat, 0 equals no repeats

# Dates to sample age distributions
datesel = c('2015-08-15', '2016-02-15')

# Combine parameters into a list
pars = list(Q=Qpars, ET=ETpars, S0=S0, C_S0=C_S0, f_thresh=f_thresh,
            spinup=spinup, datesel=datesel)

#------------------------------------------------------------------------------#
# Run model
#------------------------------------------------------------------------------#
runMod = SAS_EFs(pars, data)

#------------------------------------------------------------------------------#
# Results
#------------------------------------------------------------------------------#
# Put results in original table
data$C_Q = runMod$C_Q

# Evaluate model
evalSAS(obs = data$Cout, sim = data$C_Q) # Returns the NSE, KGE and LCCC

mean(abs(data$Cout-data$C_Q), na.rm=TRUE) # Mean residual

# Vector of simulation data where there was an observation
measC_Q = ifelse(is.na(data$Cout), NA, data$C_Q)

#------------------------------------------------------------------------------#
# Plot
#------------------------------------------------------------------------------#
# Plot of results
ggplot(data, aes(x=Date)) +
  geom_line(aes(y=Cin, col='C in')) +
  geom_point(aes(y=Cout, col='C out')) +
  geom_line(aes(y=runMod$C_Q, col='Sim C_Q')) +
  geom_point(aes(y=measC_Q, col = 'Sim C_Q during meas')) +
  labs(col='', y='C') + 
  theme(legend.position = 'bottom')

# Plot PDF
pdf_df = as_tibble(runMod$age_matr) %>% 
  mutate(Time = 1:nrow(.)) %>% 
  gather(key=Date, value=Value, -Time) %>% 
  group_by(Date) %>% 
  mutate(Value = ifelse(Time <= tail(which(Value!=0), 1)-1, Value, NA)) %>% 
  filter(!is.na(Value))

p1 = ggplot(pdf_df, aes(Time, Value, fill=Date)) +
  geom_col(width=2) +
  labs(x='Time (d)', y='Frequency (1/d)') +
  theme(legend.position='none') 

# Plot CDF
p2 = pdf_df %>% 
  group_by(Date) %>% 
  mutate(Value = cumsum(Value)) %>% 
  ggplot(., aes(Time, Value, col=Date)) +
  geom_line() + ylim(0,1) +
  labs(x='Time (d)', y='Cumulative frequency') +
  theme(legend.position = c(0,1), legend.justification = c(0,1))

# Plot together
cowplot::plot_grid(p1, p2)
