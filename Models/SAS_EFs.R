# Rewrite of SAS_EFs.m

# Implementation of the age Master Equation using a modified Euler Forward
# that accounts for the presence of 'event' water in the streamflow and
# evapotranspiration

SAS_EFs = function(pars, data) {
  # pars: list of parameters needed for the SAS_EFs
  #       Q: alpha and beta
  #       ET: alpha and beta
  #       S: S0 size of active storage
  #       f_thresh: threshold when to start sampling storage
  #       C_S0: initial concentration of tracer in storage
  # !TODO, accept other functions than just Beta
  
  # Shape parameters
  parQ  = pars$Q
  parET = pars$ET
  S0    = pars$S0
  
  # Check if there needs to be a spinup
  if(pars$spinup$spin_n > 0) {
    # Spin up the data
    source('Functions/genSpin.R')
    doSpin = genSpin(pars, data)
    data = doSpin$data
    startIndex = doSpin$startIndex
  } else {
    startIndex = 1
  }
  
  # Length of timeseries
  NN = nrow(data)
  
  # Create vectors
  S_T  = rep(0, NN) # Storage
  C_ST = rep(0, NN) # Rank storage concentration 
  C_Q  = rep(0, NN) # Tracer concentration in stream
  # !TODO: Following need to be implemented
  C_ET = rep(0, NN) # ET concentration
  C_S  = rep(0, NN) # Mean storage concentration
  
  # Age matrix that stores the age distributions at each timestep
  datesel = as.Date(pars$datesel, tz='UTC')
  age_matr = matrix(nrow=NN, ncol=length(datesel))
  colnames(age_matr) = as.character(datesel)
  
  # Initial conditions
  C_Q[1]  = pars$C_S0  # Initial conc of the stream = initial conc of storage
  C_ET[1] = pars$C_S0  # ET concentration
  C_S[1]  = pars$C_S0  # Mean storage concentration
  length_s=1                 # Rank storage vector length
  S_T[length_s]  = S0        # Initial rank storage (mm)
  C_ST[length_s] = pars$C_S0 # Mean concentration of the initial rank storage
  
  # Other
  f_thresh = pars$f_thresh
  
  # Initial SAS functions Omegas
  Omega_Q  = fSAS_beta(S_T[1:length_s]/S_T[length_s], parQ)
  Omega_ET = fSAS_beta(S_T[1:length_s]/S_T[length_s], parET)
  
  # Loop
  for(i in 1:(NN-1)) {
    # --------------------------------
    # Define domain for SAS evaluation
    # --------------------------------
    age1 = max(0, data$P[i] - data$Q[i]*Omega_Q[1] - data$ET[i]*Omega_ET[1])
    dom  = (c(0, S_T[1:length_s])+age1) / (S_T[length_s] + age1)
    
    # ----------------------------------
    # Evaluate SAS functions over domain
    # ----------------------------------
    Omega_Q  = fSAS_beta(dom, pars$Q)
    Omega_ET = fSAS_beta(dom, pars$ET)
    
    # ---------------------------------
    # Solve the master equation balance
    # ---------------------------------
    S_T[1:(length_s+1)] = pmax(0,
                               c(0,S_T[1:length_s]) + data$P[i] - 
                                 data$Q[i]*Omega_Q - data$ET[i]*Omega_ET)
    # Ensure S_T is not decreasing
    for(j in 2:(length_s+1)) {
      S_T[j] = max(S_T[j], S_T[j-1]) 
    }
    
    # -------------------------------------------
    # Update solute concentration for each parcel
    # -------------------------------------------
    # Check if there is rain, otherwise set Cin to 0
    if(data$P[i]==0) data$Cin[i]=0
    # Adjust storage
    C_ST[2:(length_s+1)] = C_ST[1:length_s] # Conservative transport of the elements
    C_ST[1] = data$Cin[i]                   # Concentration of the new input
    
    # ---------------------------
    # 4: Grow vectors if required
    # ---------------------------
    if(i==1 || S_T[length_s] < f_thresh*S_T[length_s+1]) {
      length_s = length_s + 1
      
    } else {
      # Update mean concentration of the pool
      C_ST[length_s] = max(0, 
                           (C_ST[length_s+1] * (S_T[length_s+1] - S_T[length_s]) +
                              C_ST[length_s] * (S_T[length_s]-S_T[length_s-1])) /
                             (S_T[length_s+1] - S_T[length_s-1]), na.rm=TRUE)
      
      # Merge oldest elements of S_T
      S_T[length_s] = S_T[length_s+1]
      # Merge oldest values of omega functions
      Omega_Q[length_s] = Omega_Q[length_s+1]
      Omega_Q = Omega_Q[1:length_s]
      Omega_ET[length_s] = Omega_ET[length_s+1]
      Omega_ET = Omega_ET[1:length_s]
    }
    # ----------------------------
    # Compute concentrations
    # ----------------------------
    # Q
    pQ = diff(c(0, Omega_Q))
    C_Q[i+1] = (C_ST[1:length_s] %*% pQ)[1,1]
    # ET
    pET = diff(c(0, Omega_ET))
    C_ET[i+1] = (C_ST[1:length_s] %*% pET)[1,1] 
    
    # Mean storage
    pS = diff( c(0,S_T[1:length_s]/S_T[length_s]) )
    C_S[i+1] = (C_ST[1:length_s] %*% pS)[1,1] 
    
    # Store pQ in a matrix arranged with the days in columns
    if(i>=startIndex) {
      # Check if date is requested
      # Returns NA if false, otherwise returns position to slot in data
      checkDate = match(as.Date(data$Date[i]), datesel-1)
      if(!is.na(checkDate)) {
        age_matr[1:length_s, checkDate] = pQ
      }
    }
  }
  
  # Prune output based on startIndex in case of spinup
  C_Q  = C_Q[startIndex:NN]
  C_ET = C_ET[startIndex:NN]
  C_S  = C_S[startIndex:NN]
  
  # age_matr = age_matr[startIndex:NN, startIndex:NN]
  
  # Put output into list
  out = list(data = data,
             C_Q  = C_Q,
             C_ET = C_ET,
             C_S  = C_S,
             S_T  = S_T,
             age_matr = age_matr
  )
  return(out)
  # Done
}
