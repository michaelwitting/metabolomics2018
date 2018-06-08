#load required library
library(xcms)
source("R\\msLevelMerge.R")

#set number of cores
noOfCores <- 4

## Use socket based parallel processing on Windows systems
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(noOfCores)))
} else {
  register(bpstart(SnowParam(noOfCores)))
}

##########################################################################################
# prepare data
##########################################################################################
## Read the full data at once (all MS level). Advantage: easier to link
## MS1 and MS2 data and file has to be read only once.
xrms <- readMSData("data\\PestMix1_DDA.mzML", mode = "onDisk", centroided = TRUE)

#plot BPC
bpis <- chromatogram(xrms, aggregationFun = "max")
plot(bpis)

#Fluopicolide, exact mass = 381.965430576, [M+H]+ = 382.972706
eic <- chromatogram(xrms, aggregationFun = "max", mz = c(382.96, 382.98), rt = c(430,450))
plot(eic)

##########################################################################################
# perform peak detection on MS1 level
##########################################################################################
ms1cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3,30))
ms1data <- findChromPeaks(xrms, param = ms1cwp, msLevel = 1)

#get all peaks
chromPeaks <- chromPeaks(ms1data)

#isolate Fluopicolide
chromPeak <- chromPeaks[57,]

##########################################################################################
# get corresponding MS2 spectra
##########################################################################################
#get all MS2 spectra from DDA experiment
ms2spectra <- spectra(filterMsLevel(xrms, msLevel = 2))

#filter out fitting spectra
filteredMs2spectra <- getDdaMS2Scans(chromPeak, ms2spectra)

#mark position in EIC
abline(v = unlist(lapply(filteredMs2spectra, function(x) {return(x@rt)})), col = "red")

plot(unlist(lapply(filteredMs2spectra, function(x) {return(x@rt)})), unlist(lapply(filteredMs2spectra, function(x) {return(x@precursorMz)})))


