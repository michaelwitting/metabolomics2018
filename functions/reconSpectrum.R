#' This function isolates extracted ion chromatograms for a given list of m/z ratios
#'
#' @import xcms
#' @import MSnbase
#'
#' @param msData A OnDiskExp from which EICs shall be extracted
#' @param msPeaks List of peaks for which EICs shall be extracted
#' @param rtWindow Window in seconds to enlarge EICs around peak boundaries
#' @param msLevel Indicates from which MS level EICs shall be extracted
#' @param adjustedRtime TRUE or FALSE if adjusted retention times shall be used
#' @param fileIndex index of file for which EICs shall be extracted
#'
#' @return returns a data frame containing the m/z value of the candidates and their EICs.

reconSpectrum2 <- function() {

  #required data:
  # - precursor candidate
  # - MS1 data for EIC
  # - fragment candidates
  # - MS2 data for EICs
  # - cor-Cut-off
  # - pVal-Cut-off

  #sanity checks

  #do this on all files?
  # yes --> large and long computations
  # no --> requires specific file index or file name


  #extract EICs
  ms1eics <- chromatograms()


  #align EICs

  #correlate EICs

  #filter candidates

  #construct MS2 spectrum (as Spectrum2 object)

  #return Spectrum objects as list
}

#####
# make sub functions
#####

