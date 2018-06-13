#' This function plots the candidates EIC and the fragment EICs
#'
#' @param ms1df
#' @param ms2dfCor
#' @param precursor
#' @param corValues
#' @param save
#' @param path
#' @param prefix
plotCandEICs <- function(ms1df, ms2dfCor, precursor, corValues, save = FALSE, path, prefix = "") {

  #generate filename
  if(missing(path)) {
    fileName <- paste0(prefix, "Chroms_", round(precursor[1], 0), "@", round(precursor[4], 0), ".pdf")
  } else {
    fileName <- paste0(path, prefix, "Chroms_", round(precursor[1], 0), "@", round(precursor[4], 0), ".pdf")
  }

  #get maximum value for y
  currentmax <- max(ms1df$intensity, na.rm = TRUE)

  if(max(ms2dfCor$intensityCor, na.rm = TRUE) > currentmax) {
    currentmax <- max(ms2dfCor$intensityCor, na.rm = TRUE)
  }

  sample_colors <- palette(rainbow(length(unique(ms2dfCor$mz)) + 1))
  plot(ms1df$rtime, ms1df$intensity, type = "l", col = sample_colors[1], ylim = c(0,currentmax), xlab = "retention time", ylab = "intensity")
  title(fileName)

  count <- 1

  #iterate through mz and plot
  for(i in unique(ms2dfCor$mz)) {

    count <- count + 1
    eic <- ms2dfCor[which(ms2dfCor$mz == i),]

    points(eic$rtime, eic$intensityCor, type = "l", col = sample_colors[count])

  }

  legendtext <- c("MS1", paste(round(unique(ms2dfCor$mz), 4), round(corValues$corValue, 2), sep=" ; "))

  #add legend to plot
  legend("topright", legend = legendtext, col = sample_colors, text.font = 3, lty=1, cex = 0.75)

  if(save) {
    # safe a copy of plit
    dev.copy(pdf, fileName)
    dev.off()
  }
}

#' This function isolates extracted ion chromatograms for a given list of m/z ratios
#'
#' @param ms2spectrum
#' @param precursor
#' @param save
#' @param path
#' @param prefix
plotSpectrum <- function(ms2spectrum, precursor, save = FALSE, path, prefix = "") {

  #generate filename
  if(missing(path)) {
    fileName <- paste0(prefix, "MSMS_", round(precursor[1], 0), "@", round(precursor[4], 0), ".pdf")
  } else {
    fileName <- paste0(path, prefix, "MSMS_", round(precursor[1], 0), "@", round(precursor[4], 0), ".pdf")
  }

  #create simple stick plot
  plot(ms2spectrum$mz, ms2spectrum$intensity, type = "h", xlab = "m/z", ylab = "intesity")
  title(fileName)

  if(save) {
    # safe a copy of plit
    dev.copy(pdf, fileName)
    dev.off()
  }
}
