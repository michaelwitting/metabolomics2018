#' This function isolates extracted ion chromatograms for a given list of m/z ratios
#'
#' @import xcms
#' @import MSnbase
#'
#' @param msData A OnDiskExp from which EICs shall be extracted
#' @param msPeaks List of peaks for which EICs shall be extracted
#'
#' @return returns a data frame containing the m/z value of the candidates and their EICs.
getPocketNumber <- function(swathPockets, ms1peak) {

  #get fitting SWATH window
  pocketNumber <- as.numeric(row.names(swathPockets[which(ms1peak[1] > swathPockets$start & ms1peak[1] < swathPockets$end),]))
  if (length(pocketNumber)>0 & is.numeric(pocketNumber)) {

    return(pocketNumber)

    #use higher pocketnumber for peaks on border/overlap of swath pockets
    if(length(pocketNumber) > 1) {
      pocketNumber <- pocketNumber[length(pocketNumber)]
      return(pocketNumber)
    }
  } else {
    warning("Problem!!!")
  }
}

#' This function isolates extracted ion chromatograms for a given list of m/z ratios
#'
#' @param tbl data frame with mz, rtmin and rtmax for each peak that shall be picked
#' @param msnExp OnDiskExp, needed to get scan number for rtmin and rtmax
#' @param error numeric error for isolation in m/z direction
#' @param errortype either abs or ppm (ppm not implemented yet)
#'
#' @return returns a data frame containing the m/z value of the candidates and their EICs.
tbl2ROI <- function(tbl, msnExp, error, errortype = "abs") {

  #generate empty list for ROIs
  ROIList <- list()

  #depended on error type lower and upper limit are calculated
  if(errortype == "abs") {

    #iterate through data frame with m/z rt pairs
    for(i in 1:nrow(tbl)) {

      ROIList <- c(ROIList, list(list(mz=tbl$mz[i],
                                      mzmin=tbl$mz[i] - error,
                                      mzmax=tbl$mz[i] + error,
                                      scmin=which.min(abs(tbl$rtMin[i]-MSnbase::rtime(msnExp)))[[1]],
                                      scmax=which.min(abs(tbl$rtMax[i]-MSnbase::rtime(msnExp)))[[1]],
                                      length=-1,
                                      intensity=-1)))

    }

    return(ROIList)

  } else if(errortype == "ppm") {

    #iterate through data frame with m/z rt pairs
    for(i in 1:nrow(tbl)) {

      ROIList <- c(ROIList, list(list(mz=tbl$mz[i],
                                      mzmin=tbl$mz[i] - error,
                                      mzmax=tbl$mz[i] + error,
                                      scmin=which.min(abs(tbl$rtMin[i]-MSnbase::rtime(msnExp)))[[1]],
                                      scmax=which.min(abs(tbl$rtMax[i]-MSnbase::rtime(msnExp)))[[1]],
                                      length=-1,
                                      intensity=-1)))

    }


    return(ROIList)

  } else {
    stop("wrong error type!!!")
  }
}

#' This function isolates extracted ion chromatograms for a given list of m/z ratios
#'
#' @param tbl data frame with mz, rtmin and rtmax for each peak that shall be picked
#' @param msnExp OnDiskExp, needed to get scan number for rtmin and rtmax
#' @param error numeric error for isolation in m/z direction
#' @param errortype either abs or ppm (ppm not implemented yet)
#'
#' @return returns a data frame containing the m/z value of the candidates and their EICs.
tbl2ROImzOnly <- function(tbl, msnExp, error, errortype = "abs") {

  #generate empty list for ROIs
  ROIList <- list()

  #depended on error type lower and upper limit are calculated
  if(errortype == "abs") {

    #iterate through data frame with m/z rt pairs
    for(i in 1:nrow(tbl)) {

      ROIList <- c(ROIList, list(list(mz=tbl$mz[i],
                                      mzmin=tbl$mz[i] - error,
                                      mzmax=tbl$mz[i] + error,
                                      scmin=min(MSnbase::rtime(msnExp)),
                                      scmax=max(MSnbase::rtime(msnExp)),
                                      length=-1,
                                      intensity=-1)))

    }

    return(ROIList)

  } else if(errortype == "ppm") {

    #iterate through data frame with m/z rt pairs
    for(i in 1:nrow(tbl)) {

      ROIList <- c(ROIList, list(list(mz=tbl$mz[i],
                                      mzmin=tbl$mz[i] - error,
                                      mzmax=tbl$mz[i] + error,
                                      scmin=min(MSnbase::rtime(msnExp)),
                                      scmax=max(MSnbase::rtime(msnExp)),
                                      length=-1,
                                      intensity=-1)))

    }


    return(ROIList)

  } else {
    stop("wrong error type!!!")
  }
}
