#' This function searches fragment candidates for a specified neutral loss
#'
#' @param precursor putative precursor peak
#' @param ms2peaks chromatographic peaks from SWATH window fitting precursor mass
#' @param neutralLossMass mass of neutral loss that shall be searched, e.g. 18 for water loss
#' @param mzErrorType type of error for window calculation
#' @param mzTol m/z tolerance for search window
#' @param rtTol rt tolerance for search window
#'
#' @return returns all list of candidate peaks
searchNLcandidates <- function(precursor, ms2peaks, neutralLossMass, mzErrorType = "abs", mzTol, rtTol) {

  #mass to search for in MS2 data
  searchMz = precursor[1] - neutralLossMass

  #subset data to fitting peaks
  if(mzErrorType == "abs") {

    #select candidate peaks
    ms2peaksCand <- ms2peaks[which(abs(ms2peaks[,4] - precursor[4]) < rtTol & abs(ms2peaks[,1] - searchMz) < mzTol),]

  } else if(mzErrorType == "ppm") {

    #select candidate peaks
    ms2peaksCand <- ms2peaks[which(abs(ms2peaks[,4] - precursor[4]) < rtTol & abs(ms2peaks[,1] - searchMz) / searchMz * 10^6 < mzTol),]

  } else {

    stop("Wrong error type supplied")

  }

  return(ms2peaksCand)

}

#' This function searches fragment candidates in a rt window around the specified precursor
#'
#' @param precursor putative precursor peak
#' @param ms2peaks chromatographic peaks from SWATH window fitting precursor mass
#' @param rtTol rt tolerance for search window
#'
#' @return returns all list of candidate peaks
searchMultiplePIcandidates <- function(precursor, ms2peaks, rtTol) {

  #select candidate peaks
  ms2peaksCand <- ms2peaks[which(abs(ms2peaks[,4] - precursor[4]) < rtTol),]

  #return candidates
  return(ms2peaksCand)

}

#' This function searches fragment candidates in a rt window around the specified precursor
#'
#' @param precursor putative precursor peak
#' @param ms2peaks chromatographic peaks from SWATH window fitting precursor mass
#' @param mz mz ratio of specific product ion that shall be searched
#' @param mzErrorType type of error for window calculation
#' @param mzTol m/z tolerance for search window
#' @param rtTol rt tolerance for search window
#'
#' @return returns all list of candidate peaks
searchSinglePIcandidates <- function(precursor, ms2peaks, mz, mzErrorType = "abs", mzTol, rtTol) {

  #mass to search for in MS2 data
  searchMz = mz

  #subset data to fitting peaks
  if(mzErrorType == "abs") {

    #select candidate peaks
    ms2peaksCand <- ms2peaks[which(abs(ms2peaks[,4] - rtime) < rttol & abs(ms2peaks[,1] - searchMz) < mztol),]

  } else if(mzErrorType == "ppm") {

    #select candidate peaks
    ms2peaksCand <- ms2peaks[which(abs(ms2peaks[,4] - rtime) < rttol & abs(ms2peaks[,1] - searchMz) / searchMz * 10^6 < mztol),]

  } else {

    stop("Wrong error type supplied")

  }

  return(ms2peaksCand)

}
