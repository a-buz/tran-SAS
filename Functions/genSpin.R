# Generate spinup data
genSpin = function(pars, data) {
  spinStart = pars$spinup$spinStart
  spinEnd   = pars$spinup$spinEnd
  spin_n    = pars$spinup$spin_n
  
  # Check if data is hourly or daily
  if(inherits(data$Date, "POSIXct")) {
    spinStart = as.POSIXct(spinStart, tz='UTC')
    spinEnd   = as.POSIXct(paste(spinEnd, '23:00'), tz='UTC')
  } else {
    spinStart = as.Date(spinStart, tz='UTC')
    spinEnd   = as.Date(spinEnd,   tz='UTC')
  }
  
  # Get subset to spinup
  spin     = data[data$Date >= spinStart & data$Date <= spinEnd,]
  # Repeat n amount of times
  spinning = do.call("rbind", replicate(spin_n, spin, simplify=FALSE))
  # Bind onto dataset
  spun     = rbind(spinning, data)
  
  # Find where we want to actually simulate
  startIndex = nrow(spinning)+1
  
  out = list(data=spun, startIndex=startIndex)
  return(out)
}