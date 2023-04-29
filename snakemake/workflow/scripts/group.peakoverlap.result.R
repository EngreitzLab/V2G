suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(optparse))
option.list <- list(
	make_option("--outDir", type="character"),
        make_option("--cs_peak_overlap_summary", type="character"),
        make_option("--PeaksOverlapFull_Peaks", type="character"),
        make_option("--grouped_celltype_table", type="character"),
        make_option("--cellGroup", type="character"),
        make_option("--source", type="character"),
        make_option("--PIP", type="numeric"),
	make_option("--helperFunctions", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
setwd(opt$outDir)
source(opt$helperFunctions)
cs_peak_overlap_summary <- read.table(opt$cs_peak_overlap_summary, header=TRUE, stringsAsFactors=FALSE) 
cellTags_df <- read.table(opt$grouped_celltype_table, header=TRUE, stringsAsFactors=F)
peak_info <- read.table(opt$PeaksOverlapFull_Peaks, header=TRUE, stringsAsFactors=FALSE) 
if (!any(is.na(peak_info$PIP))) {
	peak_info <- peak_info %>% filter(PosteriorProb>=opt$PIP)
}


if (opt$cellGroup %in% colnames(cs_peak_overlap_summary)){
	cs_peak_overlap_summary <- cs_peak_overlap_summary %>% select(c(CredibleSet, !!sym(opt$cellGroup)))

	if (opt$source != ""){
    		cs_peak_overlap_summary <- cs_peak_overlap_summary %>% mutate(source=opt$source) %>% separate(CredibleSet, c("remove", "LeadVariant"), sep="-", remove=FALSE) %>% select(-remove) %>% relocate(LeadVariant, .after=source)
	} else {
    		cs_peak_overlap_summary <- cs_peak_overlap_summary %>% separate(CredibleSet, c("source", "LeadVariant"), sep="-", remove=FALSE)%>% relocate(LeadVariant, .after=source)
	}

	rownames(cellTags_df) <- cellTags_df$Cellgroup
	celltype_tags <- unlist(strsplit(as.character(cellTags_df[opt$cellGroup, "Celltypes"]), split="|", fixed=TRUE))

	peak_info_EC_select <- peak_info %>% select(c("PeakChr","PeakStart","PeakEnd", "CredibleSet", "CellType","PosteriorProb")) %>% filter(greplany(celltype_tags, CellType))

	combined.credibleset.peak <- left_join(cs_peak_overlap_summary, peak_info_EC_select, by=c("CredibleSet"="CredibleSet")) %>% distinct(.keep_all=TRUE)

	write.table(combined.credibleset.peak, quote=FALSE, sep="\t", paste0("CredibleSetPeakOverlapSummary.with", opt$cellGroup, "peakinfo.tsv"))

	peak_info <- peak_info %>% select(c("PeakChr","PeakStart","PeakEnd", "CredibleSet", "CellType","PosteriorProb"))%>% mutate(isNotCelltype=ifelse(greplany(celltype_tags, CellType), 0, 1))
	celltype_specific_colnames <- paste0("OnlyIntersectsWith", opt$cellGroup, "Peak")

	peak_info <- peak_info %>% group_by(CredibleSet) %>% mutate(ContainsOtherCellTypes=ifelse(sum(isNotCelltype), TRUE, FALSE)) %>% ungroup() %>% filter(!ContainsOtherCellTypes) %>% as.data.frame()%>% select(CredibleSet, ContainsOtherCellTypes) %>% dplyr::rename(!!sym(celltype_specific_colnames):=ContainsOtherCellTypes) %>% mutate(!!sym(celltype_specific_colnames):="TRUE")

	combined.celltype.credibleset.peak <- left_join(cs_peak_overlap_summary, peak_info, by=c("CredibleSet"="CredibleSet")) %>% distinct(.keep_all=TRUE)

	write.table(combined.celltype.credibleset.peak, quote=FALSE, sep="\t", paste0("CredibleSetPeakOverlapSummary.withpeakinfo.only",opt$cellGroup,".tsv"))
} else {
	combined.credibleset.peak <- data.frame(matrix(ncol=1, nrow=0))
        write.table(combined.credibleset.peak, quote=FALSE, sep="\t", paste0("CredibleSetPeakOverlapSummary.with", opt$cellGroup, "peakinfo.tsv"))
	combined.celltype.credibleset.peak <- data.frame(matrix(ncol=1, nrow=0))
	write.table(combined.celltype.credibleset.peak, quote=FALSE, sep="\t", paste0("CredibleSetPeakOverlapSummary.withpeakinfo.only",opt$cellGroup,".tsv"))	
}

