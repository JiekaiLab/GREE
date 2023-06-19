#' Title
#'
#' @param genes A `GRanges` object that contains genes information.
#' @param binByGroup binByGroup Matrix
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from the BinByCellMatrix.
#' @param useGeneBoundaries A boolean value indicating whether gene boundaries should be employed during gene activity score
#' @param extendUpstream The minimum and maximum number of basepairs upstream of the transcription start site to consider for gene
#' activity score calculation.
#' @param extendDownstream The minimum and maximum number of basepairs downstream of the transcription termination site to consider for
#' gene activity score calculation.
#'
#' @return GeneByGroup
#' @export

getGeneByGroupMatrix <- function(genes,binByGroup,
                                 excludeChr = c("chrY","chrM"),
                                 useGeneBoundaries = TRUE,
                                 extendUpstream = 2000,
                                 extendDownstream = 2000
){
  # geneRegions <- genes[BiocGenerics::which(seqnames(genes) %bcni% excludeChr)]
  mask = !(S4Vectors::match(seqnames(genes),excludeChr,  nomatch = 0) > 0)
  geneRegions <- genes[BiocGenerics::which(mask)]

  seqlevels(geneRegions) <- as.character(unique(seqnames(geneRegions)))
  geneRegions <- geneRegions[!is.na(mcols(geneRegions)$symbol)]

  #分染色体计算gene的区域:----------------------------------------------------------------------------------
  geneRegions_list <- split(geneRegions, seqnames(geneRegions))
  geneRegions_list <- lapply(geneRegions_list, function(x){
    mcols(x)$idx <- seq_along(x)
    return(x)
  })

  extendedGeneRegion_list = lapply(seq_along(geneRegions_list), function(z){
    #Get Gene Starts
    geneRegionz <- geneRegions_list[[z]]
    geneRegionz <- geneRegionz[order(geneRegionz$idx)]
    chrz <- paste0(unique(seqnames(geneRegionz)))

    #Time to Overlap Gene Windows
    if(useGeneBoundaries){
      geneStart <- start(resize(geneRegionz, 1, "start"))
      geneEnd <- start(resize(geneRegionz, 1, "end"))

      pminGene <- pmin(geneStart, geneEnd)
      pmaxGene <- pmax(geneStart, geneEnd)

      idxMinus <- BiocGenerics::which(strand(geneRegionz) != "-")

      pForward <- rep(extendUpstream, length(pminGene))
      pForward[idxMinus] <- rep(extendDownstream, length(idxMinus))

      pReverse <- rep(extendDownstream, length(pminGene))
      pReverse[idxMinus] <- rep(extendDownstream, length(idxMinus))

      ################################################################
      #We will test when genes pass by another gene promoter
      ################################################################
      #Start of Range is based on the max observed gene ranged <- direction
      s <- pmax(
        c(1, pmaxGene[-length(pmaxGene)]),
        pminGene - pReverse
      )
      s <- pmin(pminGene, s)

      #End of Range is based on the max observed gene ranged -> direction
      e <- pmin(
        c(pminGene[-1], pmaxGene[length(pmaxGene)] + pForward[length(pmaxGene)]),
        pmaxGene + pForward
      )
      e <- pmax(pmaxGene, e)

      extendedGeneRegionz <- IRanges(start = s, end = e)

      idx1 <- which(pminGene < start(extendedGeneRegionz))
      if(length(idx1) > 0){
        stop("Error in gene boundaries minError")
      }

      idx2 <- which(pmaxGene > end(extendedGeneRegionz))
      if(length(idx2) > 0){
        stop("Error in gene boundaries maxError")
      }

    }else{
      extendedGeneRegionz <- ranges(suppressWarnings(extendGR(geneRegionz, upstream = max(extendUpstream), downstream = max(extendDownstream))))
    }

    extendedGeneRegionz = GRanges(seqnames = chrz,
                                  ranges = extendedGeneRegionz,
                                  symbol = geneRegionz$symbol)

  })


  extendedGeneRegion = extendedGeneRegion_list[[1]]
  for(i in 2:length(extendedGeneRegion_list)){
    extendedGeneRegion = suppressWarnings(c(extendedGeneRegion,extendedGeneRegion_list[[i]]))
  }

  #基于binByGroup得到bin的位置信息：
  bin = rownames(binByGroup)
  bin = stringr::str_split(bin,"_",simplify = T)
  bin_gr = GRanges(seqnames = bin[,1],
                   ranges = IRanges(start = as.numeric(bin[,2]),
                                    end = as.numeric(bin[,3])))

  tmp <- suppressWarnings(findOverlaps(extendedGeneRegion, bin_gr))
  x = rep(1,length(tmp))

  #Creating Sparse Matrix
  tmp <- Matrix::sparseMatrix(
    i = queryHits(tmp),
    j = subjectHits(tmp),
    x = x,
    dims = c(length(extendedGeneRegion), length(bin_gr))
  )

  #计算Gene Scores
  GeneByGroup <- tmp %*% as.matrix(binByGroup)
  rownames(GeneByGroup) = extendedGeneRegion$symbol

  return(GeneByGroup)

}
