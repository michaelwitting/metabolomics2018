source("functions\\eicFunctions.R")
source("functions\\fragmentCandSearch.R")
source("functions\\plottingFunctions.R")
source("functions\\reconSpectrum.R")
source("functions\\spectraFunctions.R")
source("functions\\utils.R")

## MS1 and MS2 data and file has to be read only once.
xrms <- readMSData("data\\PestMix1_SWATH.mzML", mode = "onDisk", centroided = TRUE)

#plot BPC
bpis <- chromatogram(xrms, aggregationFun = "max")
plot(bpis)

# filter out the MS1 & MS2 level
xrms_ms1 <- filterMsLevel(xrms, msLevel = 1)
xrms_ms2 <- filterMsLevel(xrms, msLevel = 2)

bpis <- chromatogram(xrms_ms1, aggregationFun = "max")
plot(bpis)

# split MS2 data into SWATH pockets
# Use the fData method to access the "feature data" (i.e. the mzML header)
precursors <- fData(xrms_ms2)$precursorMZ

#split according to different SWATH windows
xrms_ms2_swath <- split(xrms_ms2, f = as.integer(precursors))

##########################################################################################
# Peak detection in MS1 & MS2
##########################################################################################
## Set parameter for Centwave peak detection
ms1cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3,30))
ms2cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3,30))

# find peaks in MS1 level
ms1data <- findChromPeaks(xrms_ms1, param = ms1cwp, msLevel = 1)

# find peaks in MS2 level
# Peak detection on all; explicitely tell findChromPeaks to run on MS2
ms2data <- lapply(xrms_ms2_swath, FUN = findChromPeaks, param = ms2cwp, msLevel = 2)

# get all MS1 peaks
ms1chromPeaks <- chromPeaks(ms1data)

# Fluopicolide, exact mass = 381.965430576, [M+H]+ = 382.972706
eic <- chromatogram(xrms_ms1, aggregationFun = "max", mz = c(382.96, 382.98), rt = c(430,450))
plot(eic)

# isolate Fluopicolide
chromPeak <- ms1chromPeaks[41,]

##########################################################################################
# reconstruct isotope pattern + adduct, in source fragments etc
##########################################################################################
# get fragment candidates
ms1Cand <- searchMultiplePIcandidates(chromPeak, ms1chromPeaks, 5)

if(nrow(ms1Cand) > 0) {
  
  # get EICs for MS1 precursor and MS2 fragment candidates
  ms1df <- getEics(ms1data, chromPeak, rtwindow = 5, msLevel = 1)
  ms2df <- getEics(ms1data, ms1Cand, rtwindow = 5, msLevel = 1)
  
  # Align MS2 data with MS1
  ms2dfCor <- alignMsLevel(ms1df, ms2df)
  
  # correlate MS1 and MS2
  corValues <- correlateMsLevel(ms1df, ms2dfCor)
  
  # plot EICs
  plotCandEICs(ms1df, ms2dfCor, chromPeak, corValues)
  
  # reconstruct MS2 spectrum
  ms1spectrum <- reconSpectrum(chromPeak, ms2dfCor, corValues)
  plotSpectrum(ms1spectrum, chromPeak)
  
}

##########################################################################################
# reconstruct MS2 spectrum from SWATH data
##########################################################################################
# get peaks in fitting SWATH pocket
ms2chromPeaks <- chromPeaks(ms2data[[7]])

# get fragment candidates
fragCand <- searchMultiplePIcandidates(chromPeak, ms2chromPeaks, 5)

if(nrow(fragCand) > 0) {
  
  # get EICs for MS1 precursor and MS2 fragment candidates
  ms1df <- getEics(ms1data, chromPeak, rtwindow = 5, msLevel = 1)
  ms2df <- getEics(ms2data[[7]], fragCand, rtwindow = 5, msLevel = 2)
  
  # Align MS2 data with MS1
  ms2dfCor <- alignMsLevel(ms1df, ms2df)
  
  # correlate MS1 and MS2
  corValues <- correlateMsLevel(ms1df, ms2dfCor)
  
  # plot EICs
  plotCandEICs(ms1df, ms2dfCor, chromPeak, corValues)
  
  # reconstruct MS2 spectrum
  ms2spectrum <- reconSpectrum(chromPeak, ms2dfCor, corValues)
  plotSpectrum(ms2spectrum, chromPeak)
  
}
