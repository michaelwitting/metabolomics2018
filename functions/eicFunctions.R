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
getEics <- function(msData, msPeaks, rtwindow, msLevel, adjustedRtime = TRUE, fileIndex) {

  #create empty data frame for EICs
  eicDf <- data.frame()

  #get number of files
  noOfFiles <- 1
  noOfChroms <- 1

  if(!is.null(nrow(msPeaks))) {

    msPeaks[,5] <-  msPeaks[,5] - rtwindow
    msPeaks[,6] <-  msPeaks[,6] + rtwindow
    bpis <- chromatogram(msData, mz = msPeaks[,2:3], rt = msPeaks[,5:6], msLevel = msLevel)

    print(bpis)

    #depend no number of files procedure is different
    if(noOfFiles > 1) {

      #add all peaks for file
      print(noOfFiles)
      print(noOfChroms)

      for(i in 1:noOfChroms) {

        print(msPeaks[i,1])
        eicDf <- rbind.data.frame(eicDf, cbind.data.frame(mz = msPeaks[i,1],
                                                          rtime = bpis[i,fileIndex]@rtime,
                                                          intensity = bpis[i,fileIndex]@intensity))
      }

    } else {

      #add all peaks for file
      for(i in 1:length(bpis)) {
        eicDf <- rbind.data.frame(eicDf, cbind.data.frame(mz = msPeaks[i,1],
                                                          rtime = bpis[i]@rtime,
                                                          intensity = bpis[i]@intensity))
      }
    }



  } else {

    msPeaks[5] <-  msPeaks[5] - rtwindow
    msPeaks[6] <-  msPeaks[6] + rtwindow
    bpis <- chromatogram(msData, mz = msPeaks[2:3], rt = msPeaks[5:6], msLevel = msLevel)

    #depend no number of files procedure is different
    if(noOfFiles > 1) {

      #add all peaks for file
      for(i in 1:length(bpis)) {
        eicDf <- rbind.data.frame(eicDf, cbind.data.frame(mz = msPeaks[1],
                                                          rtime = bpis[1,fileIndex]@rtime,
                                                          intensity = bpis[1,fileIndex]@intensity))
      }

    } else {

      #add all peaks for file
      for(i in 1:length(bpis)) {
        eicDf <- rbind.data.frame(eicDf, cbind.data.frame(mz = msPeaks[1],
                                                          rtime = bpis[i]@rtime,
                                                          intensity = bpis[i]@intensity))
      }
    }

  }

  # replace NAs with 0 in EICs
  eicDf[is.na(eicDf)] <- 0

  #return data frame
  return(eicDf)

}

#' This function aligns the MS1 and MS2 traces to overcome the lag between them based on the different dwell times
#'
#' @param ms1df Data frame with the EIC for the precursor
#' @param ms2df data frome with the EICs for all fragment candidates
#'
#' @return returns a data frame containing the m/z value of the candidates and their EICs corrected according to the MS1 EICs
alignMsLevel <- function(ms1df, ms2df) {

  #create data frame to store corrected EICs
  ms2dfCor <- data.frame()

  #iterate through all MS2 candidate EICs
  for(i in unique(ms2df$mz)) {

    #isolate individual EIC to work with
    candEIC <- ms2df[which(ms2df$mz == i),]

    #use linear approximation to get values at RTs of MS1
    candEICCor <- approx(candEIC$rtime, candEIC$intensity, ms1df$rtime)
    rtimeCor <- candEICCor$x
    intensityCor <- candEICCor$y

    #add to result data frame
    ms2dfCor <- rbind.data.frame(ms2dfCor, cbind.data.frame(mz = i,
                                                            rtimeCor = rtimeCor,
                                                            intensityCor = intensityCor))
  }

  # replace NAs with 0 in EICs
  ms2dfCor[is.na(ms2dfCor)] <- 0

  #return new data frame
  return(ms2dfCor)

}

#' This function correlates the MS1 and MS2 traces and returns a data frame with the pearson correlation coefficients
#'
#' @param ms1df Data frame with the EIC for the precursor
#' @param ms2dfCor data frome with the EICs for all fragment candidates
#'
#' @return a data frame containing m/z ratios and correlation coefficients
correlateMsLevel <- function(ms1df, ms2dfCor) {

  #create data frame for correlation values
  corDf <- data.frame()

  #iterate through all MS2 candidate EICs
  for(i in unique(ms2dfCor$mz)) {

    #isolate individual corrected EIC to work with
    corEIC <- ms2dfCor[which(ms2dfCor$mz == i),]

    #perform correlation with Ms1 level
    corValue <- cor.test(ms1df$intensity, corEIC$intensityCor)

    #add results to data frame
    corDf <- rbind.data.frame(corDf, cbind.data.frame(mz = i,
                                                      corValue = corValue$estimate,
                                                      pValue = corValue$p.value))

  }

  #return new data frame
  return(corDf)

}
