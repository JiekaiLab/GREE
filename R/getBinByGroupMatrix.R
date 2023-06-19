#' Get the BinByGroup matrix
#'
#' @param binByCell BinByCellMatrix
#' @param cell_anno Cell meta with cell group information
#' @param n_sample Number of cells sampled in each group
#' @param group_col The column name of the cell group in cell_anno
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#'
#' @return binByGroup_norm
#' @export

getBinByGroupMatrix <- function(binByCell,cell_anno,
                                group_col="group",
                                n_sample=1000,
                                seed=1000
){
  #对于大于n_sample的group随机抽取n_sample个，小于的升采样到n_sample
  cellNames <- rownames(cell_anno)
  group = names(table(cell_anno[,group_col]))
  group_cellNum = table(cell_anno[,group_col])

  set.seed(seed)
  cell_list = list()
  binByGroup = data.frame(row.names = rownames(binByCell))
  for(i in names(group_cellNum)){
    cell = cellNames[cell_anno[,group_col]==i]
    if(group_cellNum[i]>n_sample){
      cell = sample(cell,n_sample)
    }else{
      times = trunc(n_sample/length(cell))
      cell = c(rep(cell,times),sample(cell,n_sample-times*length(cell)))
    }

    binByGroup[,i] = BiocGenerics::rowSums(binByCell[,cell])
    # mask = match(cell,colnames(binByCell))
    # binByGroup[,i] = BiocGenerics::rowSums(binByCell[,mask])

    cell_list[[i]] = cell
  }

  #只考虑筛选后的片段数来进行标准化：
  groupNormFactors <- colSums(binByGroup)
  scaleFactors <- 10^7 / groupNormFactors #Scale with Norm Factors

  binByGroup_norm = binByGroup
  for(i in group){
    binByGroup_norm[,i] = binByGroup[,i]*scaleFactors[i]
  }

  return(binByGroup_norm)
}
