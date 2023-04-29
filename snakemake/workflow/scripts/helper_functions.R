## Jesse Engreitz
## August 20, 2012
## R utilities
require(GenomicRanges)
sortChromosomeNames <- function(chrs) {
  x <- gsub("chr","",chrs)
  suppressWarnings(result <- c(paste("chr", sort(na.omit(as.numeric(x))), sep=''),
              paste("chr", sort(x[is.na(as.numeric(x))]), sep='')))
  return(result)
}
bed.extra.colnames <- c("name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
readBed <- function(file, extra.colnames=bed.extra.colnames, chr=NULL, sort=F, skip.chr.sorting=FALSE) {
  first.line <- readLines(file, n=1)
  if (regexpr("track", first.line) != -1) skip=1 else skip=0
  result <- read.delim(file, stringsAsFactors=T, header=F, comment.char="#", skip=skip)
  colnames(result) <- c("chr","start","end", extra.colnames)[1:ncol(result)]
  if (!is.null(chr)) result <- result[result$chr == chr,]
  if (!skip.chr.sorting) result$chr <- ordered(result$chr, levels=sortChromosomeNames(levels(result$chr)))
  if (sort) result <- result[order(result$chr, result$start, result$end),]
  return(result)
}


write.tab <- function(d, file, row.names=F, col.names=T, sep='\t', quote=F, ...) {
  write.table(d, file, row.names=row.names, col.names=col.names, sep=sep, quote=quote, ...)
}


expandBed <- function(bed.df, upstream, downstream) {
  if (!("strand" %in% colnames(bed.df)) | !all(bed.df$strand %in% c("-","+"))) {
    cat("Processing as unstranded\n")
    bed.df$start <- bed.df$start - upstream
    bed.df$end <- bed.df$end + downstream
  } else {
    cat("Processing as stranded\n")
    ## STRANDED
    pos <- bed.df$strand == "+"
    neg <- !pos

    bed.df$start[pos] <- bed.df$start[pos] - upstream
    bed.df$end[pos] <- bed.df$end[pos] + downstream

    bed.df$start[neg] <- bed.df$start[neg] - downstream
    bed.df$end[neg] <- bed.df$end[neg] + upstream
  }
  return(bed.df)
}

GRangesFromBed <- function(bed, ...) {
  #cat("Note:  GRangesFromBed ADJUSTS start coordinate by 1 to match 0-based BED coordinates, unless start == end (different from behavior prior to 7/8/17)\n")
  bed$start <- bed$start + 1
  bed$start[bed$start > bed$end] <- bed$end[bed$start > bed$end]
  g <- GRanges(bed[,1], IRanges(bed[,2], bed[,3]), ...)
  if (ncol(bed) > 3) names(g) <- bed[,4]
  if (ncol(bed) > 4) mcols(g) <- bed[,5]
  if (ncol(bed) > 5) {
    bed[,6] <- as.character(as.matrix(bed[,6]))
    bed[bed[,6] == ".",6] <- "+"   ## GRanges doesn't like "." as strand
    strand(g) <- bed[,6]
  }
  if (ncol(bed) > 6) mcols(g) <- bed[,c(4,7:ncol(bed))]
  return(g)
}

annotateVariants <- function(df, variant.gr, promoters.gr, include.names=F) {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(VariantAnnotation)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
  library(org.Hs.eg.db)
  codingVariants <- locateVariants(variant.gr, txdb, CodingVariants())
  df$Coding <- 1:nrow(df) %in% codingVariants$QUERYID  ## Note that this only means it overlaps CDS, not that it's missense
  if (include.names){
    df$CodingVariantGene <- NA
    keys <- unique(na.omit(codingVariants$GENEID))
    if (length(keys) > 0) {
      symbols <- select(org.Hs.eg.db, keys=keys, keytype="ENTREZID",columns="SYMBOL")
      df$CodingVariantGene[codingVariants$QUERYID] <- sapply(codingVariants$GENEID, function(id) {
        index <- which(symbols$ENTREZID == id)
        if (length(index) == 1) return(symbols$SYMBOL[index])
        else return(NA)
    })
    }
  }

  spliceVariants <- locateVariants(variant.gr, txdb, SpliceSiteVariants())
  df$SpliceSite <- 1:nrow(df) %in% spliceVariants$QUERYID
  if (include.names) {
    df$SpliceSiteVariantGene <- NA
    keys <- unique(na.omit(spliceVariants$GENEID))
    if (length(keys) > 0) {
      symbols <- select(org.Hs.eg.db, keys=keys, keytype="ENTREZID", columns="SYMBOL")
      df$SpliceSiteVariantGene[spliceVariants$QUERYID] <- sapply(spliceVariants$GENEID, function(id) {
        index <- which(symbols$ENTREZID == id)
        if (length(index) == 1) return(symbols$SYMBOL[index])
        else return(NA)
      })
    }
  }

  ixn <- findOverlaps(variant.gr, promoters.gr)
  print(ixn)
  df$Promoter <- 1:nrow(df) %in% queryHits(ixn)

  if (include.names) {
    df$PromoterVariantGene <- NA
    if (length(ixn) > 0) {
      for (i in 1:length(ixn)) {
        gene <- promoters.gr[subjectHits(ixn)[i],]$symbol
        print(gene)
        if (is.na(df$PromoterVariantGene[queryHits(ixn)[i]]))
          df$PromoterVariantGene[queryHits(ixn)[i]] <- gene
        else
          df$PromoterVariantGene[queryHits(ixn)[i]] <- paste0(unique(c(gene, strsplit(df$PromoterVariantGene[queryHits(ixn)[i]],";"))), collapse=';')
      }
    }
  }

  return(df)
}


makeCredibleSets <- function(df, genes, Source=NA, include.names=F) {
  if ("Partition" %in% colnames(df)) {
    res <- do.call(rbind, tapply(1:nrow(df), df$CredibleSet, function(i) {
       with(df[i,], data.frame(chr=as.character(as.matrix(chr[1])), start=min(position), end=max(position),CredibleSet=CredibleSet[1], nSNP=length(position), AnyCoding=any(Partition == "CDS"), AnyPromoter=any(Partition == "TSS-100bp"), AnySpliceSite=any(Partition == "SpliceSite"), Trait=Trait[1], Source=Source, BestSNP=variant[which.max(PosteriorProb)[1]]))
    }))
  } else {
    res <- do.call(rbind, tapply(1:nrow(df), df$CredibleSet, function(i) {
     with(df[i,], data.frame(chr=as.character(as.matrix(chr[1])), start=min(position), end=max(position),CredibleSet=CredibleSet[1], nSNP=length(position), AnyCoding=any(Coding), AnyPromoter=any(Promoter), AnySpliceSite=any(SpliceSite), Trait=Trait[1], Source=Source, LocusID=LocusID[1], BestSNP=variant[which.max(PosteriorProb)[1]]))
    }))
  }
  res <- merge(res, unique(data.frame(BestSNP=df$variant, BestSNPPos=df$position)), all.x=TRUE)
  res <- res[,c(2:ncol(res),1)]

  coding.genes <- subset(genes, !grepl("NR_", name))
  ixn <- as.data.frame(nearest(GRangesFromBed(with(res, data.frame(chr=chr, start=BestSNPPos, end=BestSNPPos))), GRangesFromBed(coding.genes), select="all", ignore.strand=TRUE))
  res$BestSNPNearestGene <- sapply(1:nrow(res), function(i) paste(unique(coding.genes$symbol[subset(ixn, queryHits == i)$subjectHits]), collapse=","))

  if (include.names) {
    groupGeneList <- function(strs) {
      geneList <- na.omit(strs)
      if (length(geneList) == 0) geneList <- "NA"
      uniqGenes <- geneList %>% strsplit(";") %>% unlist() %>% unique()
      paste0(uniqGenes, collapse=';')
    }

    geneLists <- df %>% group_by(CredibleSet) %>%
             summarize(CodingVariantGene=groupGeneList(CodingVariantGene),
                       SpliceSiteVariantGene=groupGeneList(SpliceSiteVariantGene),
                       PromoterVariantGene=groupGeneList(PromoterVariantGene)) %>%
             as.data.frame()

    res <- merge(res, geneLists, all.x=TRUE)
  }

  return(res)
}

greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}
