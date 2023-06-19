#' Get the BinByCell matrix
#'
#' @param frag A `GRanges` object that contains fragment information.
#' @param genes A `GRanges` object that contains genes information.
#' @param CpG A `GRanges` object that contains CpG information.
#' @param exon A `GRanges` object that contains exon information.
#' @param Blacklist A `GRanges` object that contains Blacklist information.
#' @param tileSize The size of the tiles used for binning counts in the BinByCellMatrix.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from the BinByCellMatrix.
#'
#' @return BinByCellMatrix
#' @export
#'
getBinByCellMatrix <- function(frag,genes,CpG,exon,Blacklist,
                               tileSize = 500,
                               excludeChr = c("chrY","chrM")

){
  #提取保存到"RG"列中的cell name信息：
  cellNames = unique(S4Vectors::mcols(frag)$RG)

  #去掉不想要的染色体位置，以及没有基因注释的基因组位置。
  # geneRegions <- genes[BiocGenerics::which(seqnames(genes) %bcni% excludeChr)]
  mask = !(S4Vectors::match(seqnames(genes),excludeChr,  nomatch = 0) > 0)
  geneRegions <- genes[BiocGenerics::which(mask)]
  seqlevels(geneRegions) <- as.character(unique(seqnames(geneRegions)))
  geneRegions <- geneRegions[!is.na(S4Vectors::mcols(geneRegions)$symbol)]

  #分染色体计算:----------------------------------------------------------------------------------
  geneRegions_list <- split(geneRegions, seqnames(geneRegions))
  geneRegions_list <- lapply(geneRegions_list, function(x){
    S4Vectors::mcols(x)$idx <- seq_along(x)
    return(x)
  })

  binByCell_list = lapply(seq_along(geneRegions_list), function(z){

    #Get Gene Starts
    geneRegionz <- geneRegions_list[[z]]
    geneRegionz <- geneRegionz[order(geneRegionz$idx)]
    chrz <- paste0(unique(seqnames(geneRegionz)))

    fragz = frag[seqnames(frag)==chrz]
    seqlevels(fragz) <- as.character(unique(seqnames(fragz)))

    #Read in Fragments
    fragSt <- trunc(start(fragz)/tileSize) * tileSize
    fragEd <- trunc(end(fragz)/tileSize) * tileSize

    fragBC <- rep(S4Vectors::match(S4Vectors::mcols(fragz)$RG, cellNames), 2)
    rm(fragz)
    gc()

    #Unique Inserts
    uniqIns <- sort(unique(c(fragSt,fragEd)))

    #Construct tile by cell mat!
    binByCell <- Matrix::sparseMatrix(
      i = match(c(fragSt, fragEd), uniqIns),
      j = as.vector(fragBC),
      x = rep(1,  2*length(fragSt)),
      dims = c(length(uniqIns), length(cellNames))
    )

    #Unique Tiles
    uniqueTiles <- IRanges(start = uniqIns, width = tileSize)
    #Clean Memory
    rm(uniqIns, fragSt, fragEd, fragBC)
    gc()

    rownames(binByCell) =  paste0(chrz,"_",start(uniqueTiles),"_",end(uniqueTiles))
    colnames(binByCell) <- cellNames
    return(binByCell)
  }
  )

  binByCell = binByCell_list[[1]]
  for(i in 2:length(binByCell_list)){
    binByCell = rbind(binByCell,binByCell_list[[i]])
  }

  # region_remove -----------------------------------------------------------
  # CpG = CpG[seqnames(CpG) %in% as.character(unique(seqnames(geneRegions)))]
  mask = (S4Vectors::match(seqnames(CpG),as.character(unique(seqnames(geneRegions))),  nomatch = 0) > 0)
  CpG <- CpG[BiocGenerics::which(mask)]
  seqlevels(CpG) <- as.character(unique(seqnames(CpG)))

  # exon = exon[seqnames(exon) %in% as.character(unique(seqnames(geneRegions)))]
  mask = (S4Vectors::match(seqnames(exon),as.character(unique(seqnames(geneRegions))),  nomatch = 0) > 0)
  exon <- exon[BiocGenerics::which(mask)]

  # Blacklist = Blacklist[seqnames(Blacklist) %in% as.character(unique(seqnames(geneRegions)))]
  mask = (S4Vectors::match(seqnames(Blacklist),as.character(unique(seqnames(geneRegions))),  nomatch = 0) > 0)
  Blacklist <- Blacklist[BiocGenerics::which(mask)]
  seqlevels(Blacklist) <- as.character(unique(seqnames(Blacklist)))


  extendGR <-  function(gr = NULL, upstream = NULL, downstream = NULL){
    #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
    isMinus <- BiocGenerics::which(strand(gr) == "-")
    isOther <- BiocGenerics::which(strand(gr) != "-")
    #Forward
    start(gr)[isOther] <- start(gr)[isOther] - upstream
    end(gr)[isOther] <- end(gr)[isOther] + downstream
    #Reverse
    end(gr)[isMinus] <- end(gr)[isMinus] + upstream
    start(gr)[isMinus] <- start(gr)[isMinus] - downstream
    return(gr)
  }

  TSSRegions <- resize(geneRegions, 1, "start")
  TSSRegions <- extendGR(gr = TSSRegions, upstream = 50, downstream = 100) #TSSUpstream=50,TSSDownstream=100
  o = findOverlaps(CpG,TSSRegions)
  CpG_tss = CpG[unique(queryHits(o))]

  gr_remove = c(CpG_tss,exon,Blacklist)
  gr_remove = unique(gr_remove)

  bin = stringr::str_split(rownames(binByCell),"_",simplify = T)
  bin = GRanges(seqnames = bin[,1],
                ranges = IRanges(start = as.numeric(bin[,2]),end = as.numeric(bin[,3]))
  )
  mask = !overlapsAny(bin, gr_remove)

  binByCell = binByCell[mask,]
  return(binByCell)
}

