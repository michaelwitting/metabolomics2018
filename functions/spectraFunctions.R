#' This function reconstructs a Ms2 spectrum for all fragments with a correlation value bigger than a set cut-off
#'
#' @param precursor putative precursor peak
#' @param ms2dfCor corrected MS2 EICs
#' @param corValues calculated pearson correlations
#' @param CorCutOff cut-off for pearson correlation, default = 0.9
#'
#' @return returns a data frame containing the reconstructed MS2 spectrum
reconSpectrum <- function(precursor, ms2dfCor, corValues, CorCutOff = 0.9, pValCutOff = 0.05) {

  #create empty data frame for MS2 spectrum
  spectraDf <- data.frame()

  #iterate through corrlation values and check if bigger than cut off
  for(i in corValues$mz[which(corValues$corValue > CorCutOff & corValues$pValue < pValCutOff)]) {

    #get fitting EICs
    rTime <- ms2dfCor$rtimeCor[which(ms2dfCor$mz == i)]
    intensity <- ms2dfCor$intensityCor[which(ms2dfCor$mz == i)]

    #get intensity at apex of MS1 peak
    int <- approx(rTime, intensity, precursor[4])$y

    #add to data frame
    spectraDf <- rbind.data.frame(spectraDf,
                                  cbind.data.frame(mz = i,
                                                   intensity = int,
                                                   corValue = corValues$corValue[which(corValues$mz == i)],
                                                   pValue = corValues$pValue[which(corValues$mz == i)]))

  }

  #return spectrum as data frame
  return(spectraDf)

}

#' This function creates an MS2 spectrum object from a precursor and a MS2 data frame
#'
#' @import MSnbase
#'
#' @param precursor putative precursor peak
#' @param ms2dfCor corrected MS2 EICs
#' @param corValues calculated pearson correlations
#' @param CorCutOff cut-off for pearson correlation, default = 0.9
#'
#' @return returns a data frame containing the reconstructed MS2 spectrum
makeSpectrum2 <- function() {



}
