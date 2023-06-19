#' Get the High variability bin matrix
#'
#' @param binByGroup_norm binByGroup_norm matrix
#' @param bin_max_thr Used to sieve out very low peaks, i.e. peaks with small maximum values, which may be noise.
#' @param bin_fc_thr Used to screen the bin for high variance
#'
#' @return binByGroup_norm_HVbin
#' @export

getHVbinMatrix <- function(binByGroup_norm,
                           bin_max_thr = 7, #用于筛掉很低的峰，即最大值也很小的峰，这些峰可能为噪音。
                           bin_fc_thr = 70 #用于筛选高变异的bin
){
  bin_max = matrixStats::rowMaxs(as.matrix(binByGroup_norm))
  #   bin_max_thr = quantile(bin_max,bin_max_probs)
  binByGroup_norm_HVbin = binByGroup_norm[bin_max>bin_max_thr,]

  binByGroup_norm_HVbin[binByGroup_norm_HVbin==0]=0.1

  fc = matrixStats::rowMaxs(as.matrix(binByGroup_norm_HVbin))/matrixStats::rowMins(as.matrix(binByGroup_norm_HVbin))

  fc_thr = bin_fc_thr

  binByGroup_norm_HVbin = binByGroup_norm_HVbin[fc>fc_thr,]

  return(binByGroup_norm_HVbin)
}
