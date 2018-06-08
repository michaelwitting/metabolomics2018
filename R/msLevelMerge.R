#function for DDA data
#MS2 spectra need to be list of spectrum2 objects
getDdaMS2Scans <- function(chromPeak, ms2spectra, mzTol = 0.01, rtTol = 20) {
  
  #make a new list
  filteredMs2spectra <- list()
  
  #isolate only the ones within range
  for(i in 1:length(ms2spectra)) {
    
    #isolate Ms2 spectrum
    ms2spectrum <- ms2spectra[[i]]
    
    #check if within range of peak
    if(abs(chromPeak[4] - ms2spectrum@rt) < rtTol & abs(chromPeak[1] - ms2spectrum@precursorMz) < mzTol) {
      filteredMs2spectra <- c(filteredMs2spectra, ms2spectrum)
    }
  }
  
  #return list with filtered ms2 spectra
  return(filteredMs2spectra)
}
