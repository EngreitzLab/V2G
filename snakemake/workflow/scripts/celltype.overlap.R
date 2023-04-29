suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

option.list <- list(
	make_option("--ABCOverlap", type="character", help="Path to ABCOverlap_file"),
	make_option("--PeakOverlap", type="character", help="Path to PeakOverlap_file"),
	make_option("--removeNonCoding", type="logical", help="Whether to remove non-coding genes"),
	make_option("--variantList", type="character", help="Absolute path to the list of variants"),
	make_option("--csList", type="character", help="Absolute path to the list of credible sets"), 
	make_option("--outDir", type="character", help="Path to output directory"),
	make_option("--grouped_celltype_table", type="character", help="Cell type group table"),
	make_option("--PIP", type="numeric",  help="The PIP threshold to filter the variants"),
	make_option("--helperFunctions", type="character")
	)

opt <- parse_args(OptionParser(option_list=option.list))
source(opt$helperFunctions)
setwd(opt$outDir)
abc <- read.delim(opt$ABCOverlap)
if (opt$removeNonCoding){
	abc <- abc %>% filter(!greplany(c("^LINC","-AS","^MIR","RNU","^LOC"),TargetGene))
}else{
	abc <- abc}

variants <- read.delim(opt$variantList)
all.cs <- read.delim(opt$csList)

geneRanks <- abc %>% group_by(CredibleSet,TargetGene) %>% summarise(MaxABC=max(ABC.Score)) %>% as.data.frame()
geneRanks <- geneRanks %>% group_by(CredibleSet) %>% mutate(GeneRank=dense_rank(-MaxABC)) %>% as.data.frame()
abc.ranked <- merge(abc, geneRanks)
write.table(abc.ranked, file="ABCOverlapFull.ranked.tsv",row.names=F, col.names=T, sep='\t', quote=F)

cellTags_df <- read.table(opt$grouped_celltype_table, header=TRUE, stringsAsFactors=F)
cellTags <- vector(mode="list",length=dim(cellTags_df)[1])
for (row in 1:nrow(cellTags_df)){
	entry=unlist(strsplit(as.character(cellTags_df[row, "Celltypes"]), split="|",fixed = TRUE))
	group=unlist(as.character(cellTags_df[row, "Cellgroup"]))
	cellTags[[group]] <- entry
	
} 

tmp <- abc.ranked %>% group_by(variant,CredibleSet) %>% 
  summarise(
  	AllCellTypes=paste0(CellType, collapse=','), 
  	TargetGenes=paste0(unique(TargetGene), collapse=','), 
  	TopGene=unique(TargetGene[which.max(ABC.Score)])) %>% 
  merge(all.cs) %>% arrange(CredibleSet) %>% unique() %>%  as.data.frame()

write.table(tmp, file="ABCVariantOverlapSummary.tsv",row.names=F, col.names=T, sep='\t', quote=F)


topGenes <- lapply(names(cellTags), function(category) {
	tags <- cellTags[[category]]
	if (any(is.na(abc.ranked$PosteriorProb))) {
	abc.ranked %>% filter(GeneRank == 1 & greplany(tags, CellType)) %>% select(TargetGene) %>% unique() %>% mutate(CellCategory=category)
} else {
	abc.ranked %>% filter(GeneRank == 1 & PosteriorProb >= opt$PIP & greplany(tags, CellType)) %>% select(TargetGene) %>% unique() %>% mutate(CellCategory=category)
	}
})
names(topGenes) <- names(cellTags)

topGeneTable <- do.call(rbind, topGenes) %>% spread(CellCategory,CellCategory,fill='NA')
topGeneTable[,-1:-1][topGeneTable[,-1:-1] != "NA"] <- TRUE
topGeneTable[,-1:-1][topGeneTable[,-1:-1] == "NA"] <- FALSE
write.table(topGeneTable, file="ABCTopGeneTableWithCellCategories.tsv", row.names=F, col.names=T, sep='\t', quote=F)


peakOverlap <- read.delim(opt$PeakOverlap)
csPeakOverlap <- lapply(names(cellTags), function(category) {
	tags <- cellTags[[category]]
	if (any(is.na(peakOverlap$PosteriorProb))) {
	peakOverlap %>% filter(greplany(tags, CellType)) %>% select(CredibleSet) %>% unique() %>% mutate(CellCategory=category)
} else {	
	peakOverlap %>% filter(PosteriorProb >= opt$PIP & greplany(tags, CellType)) %>% select(CredibleSet) %>% unique() %>% mutate(CellCategory=category)
}
})
names(topGenes) <- names(cellTags)
csPeakOverlap <- do.call(rbind, csPeakOverlap) %>% spread(CellCategory,CellCategory,fill='NA')
csPeakOverlap[,-1][csPeakOverlap[,-1] != "NA"] <- TRUE
csPeakOverlap[,-1][csPeakOverlap[,-1] == "NA"] <- FALSE
csPeakOverlap <- merge(csPeakOverlap, all.cs, by="CredibleSet")
write.table(csPeakOverlap, file=paste0("CredibleSetPeakOverlapSummary.tsv"), row.names=F, col.names=T, sep='\t', quote=F)
