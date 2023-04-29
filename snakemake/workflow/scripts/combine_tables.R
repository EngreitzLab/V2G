suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(optparse))

option.list <- list(
	make_option("--CSbase_df", type="character"),
	make_option("--ABCPeak_overlapping_table", type="character"),
	make_option("--UBQ_genes", type="character"),
	make_option("--LipidBloodAssociationTable", type="character"), 
	make_option("--Cellgroups", type="character"),
	make_option("--Source", type="character"),
	make_option("--OutDir", type="character"),
	make_option("--trait", type="character")
	)

opt <- parse_args(OptionParser(option_list=option.list))

# all the "MESDC1" are replaced with "TLNRD1"
CS.base.df <- read.table(opt$CSbase_df, sep = '\t', header=TRUE,stringsAsFactors=FALSE)

peak_ABC_overlapping_table <- read.table(opt$ABCPeak_overlapping_table, header=TRUE, sep = '\t', stringsAsFactors=FALSE) %>% column_to_rownames(., var = "Cellgroup")

UBQ.genes <- read.table(opt$UBQ_genes, header=FALSE,stringsAsFactors=FALSE)
lipid.blood.association <- read.table(opt$LipidBloodAssociationTable, sep =',', header=TRUE, stringsAsFactors=FALSE)
cellgroups=as.list(strsplit(opt$Cellgroups, ",")[[1]])
Source=opt$Source

#############################################################
#lipid level association
lipid.keys <- c("LDL-direct", "Triglycerides", "Cholesterol", "HDL-cholesterol", "Apolipoprotein_B", "Apolipoprotein B", "Apolipoprotein A", "Apolipoprotein_A", "HDLC", "LDLC")
#blood pressure associations
blood.keys <- c("DBP", "SBP")
#############################################################
CS.base.df <- CS.base.df %>% mutate(gene=ifelse(gene=="MESDC1","TLNRD1",gene))
#############################################################
# Finish adding additional genes 
#############################################################
#Add UbiquitousGene info
colnames(UBQ.genes) <- "UbiquitouslyExpressedGenes"
UBQ.genes$IsThisGeneUbiquitouslyExpressed <- "TRUE"
CS.base.df <- left_join(CS.base.df,UBQ.genes, by=c("gene" = "UbiquitouslyExpressedGenes"))
CS.base.df <- CS.base.df %>% mutate(IsThisGeneUbiquitouslyExpressed=ifelse(is.na(IsThisGeneUbiquitouslyExpressed),"FALSE","TRUE"))

#############################################################
# Combine with lipid level associations
lipid.blood.association.core <- lipid.blood.association[c("rsID","UKBB.PheWAS.Continuous.traits")] %>% rename("AssociatedTraits"="UKBB.PheWAS.Continuous.traits")
lipid.snps <- c()
blood.snps <- c()
for (row in 1:nrow(lipid.blood.association.core )){
	traits <- strsplit(as.character(lipid.blood.association.core[row,"AssociatedTraits"]),split="|",fixed=TRUE)[[1]]
	snp <- as.character(lipid.blood.association.core[row,"rsID"])
	test.length <- length(traits[which(traits %in% lipid.keys)])
	if (length(traits[which(traits %in% lipid.keys)]) > 0){
		lipid.snps <- c(lipid.snps, snp)
	}
	if (length(traits[which(traits %in% blood.keys)])> 0){
		blood.snps <- c(blood.snps, snp)
	}

}
CS.base.df <- CS.base.df %>% mutate(LipidLevelsAssociated=ifelse(LeadVariant %in% lipid.snps, "TRUE", "FALSE")) %>% mutate(BloodPressureAssociated=ifelse(LeadVariant %in% blood.snps, "TRUE", "FALSE"))
#############################################################
counter=0
for (cellgroup in cellgroups){
	counter = counter+1 
	#Peaks
	peaks_df <- read.table(peak_ABC_overlapping_table[cellgroup,"Peaks"],header=TRUE, sep = '\t', stringsAsFactors=FALSE)
	only_intersects_with_cellgroup_peaks_df <- read.table(peak_ABC_overlapping_table[cellgroup,"OnlyPeaks"],header=TRUE, sep = '\t', stringsAsFactors=FALSE)
	#Remove duplications caused by multiple variants in a credible sets.
	if (("CellType" %in% colnames(peaks_df)) && (paste0("OnlyIntersectsWith", cellgroup, "Peak") %in% colnames(only_intersects_with_cellgroup_peaks_df))){
		peaks.core <- peaks_df[c(cellgroup,"LeadVariant","PeakChr", "PeakStart","PeakEnd", "CellType")] %>% distinct() %>% rename(!!paste0("CellTypeOftheClosest",cellgroup,"Peak"):=CellType) %>% rename(!!paste0("IntersectsWith",cellgroup,"Peak") := !!cellgroup) 
		CS.base.df.peaks <- left_join(CS.base.df, peaks.core, by=c("LeadVariant"="LeadVariant"))
		CS.base.df.peaks <- CS.base.df.peaks %>% rowwise() %>% mutate(!!sym(paste0(cellgroup, "_PeakCenter")):=median(c(PeakStart,PeakEnd))) 
		CS.base.df.peaks <- CS.base.df.peaks %>% rowwise() %>% mutate(!!sym(paste0("DistanceToNearest",cellgroup, "PeakWithVariant")):=!!sym(paste0(cellgroup, "_PeakCenter"))-as.numeric(TSS_position)) %>% mutate(!!sym(paste0("absDistanceToNearest",cellgroup,"PeakWithVariant")):=abs(!!sym(paste0("DistanceToNearest", cellgroup, "PeakWithVariant"))))
		# Pick the gene with the shortest distance to nearest gene with cellgroup peak
		# get the gene with the minimum distance
		CS.base.df.peaks <- CS.base.df.peaks %>% group_by(Source,LeadVariant, gene) %>% filter(!!sym(paste0("absDistanceToNearest",cellgroup,"PeakWithVariant")) == min(!!sym(paste0("absDistanceToNearest",cellgroup,"PeakWithVariant")), na.rm=TRUE) |all(is.na(!!sym(paste0("absDistanceToNearest",cellgroup,"PeakWithVariant"))))) %>% ungroup()
		CS.base.df.peaks <- CS.base.df.peaks %>% group_by(Source,LeadVariant) %>% mutate(!!sym(paste0("RankOfDistanceTo", cellgroup, "PeakWithVariant")):=dense_rank(!!sym(paste0("absDistanceToNearest",cellgroup,"PeakWithVariant")))) %>% ungroup()

		only_intersects_with_cellgroup_peaks.core <- only_intersects_with_cellgroup_peaks_df %>% select(c(LeadVariant,!!sym(paste0("OnlyIntersectsWith", cellgroup,"Peak")))) %>% filter(!is.na(!!sym(paste0("OnlyIntersectsWith", cellgroup,"Peak"))))

		CS.base.df <- left_join(CS.base.df.peaks,only_intersects_with_cellgroup_peaks.core, by=c("LeadVariant"="LeadVariant")) %>% relocate(!!sym(paste0(cellgroup, "_PeakCenter")), .after=!!sym(paste0("CellTypeOftheClosest",cellgroup,"Peak"))) %>% relocate(!!sym(paste0("OnlyIntersectsWith", cellgroup,"Peak")),.after=!!sym(paste0("IntersectsWith",cellgroup,"Peak"))) 
} 
	#ABC
	ABC_df <- read.table(peak_ABC_overlapping_table[cellgroup,"ABC"],header=TRUE, sep = '\t', stringsAsFactors=FALSE)
	if (nrow(ABC_df)!=0){
		ABC.df.core <- ABC_df %>% select(c(TargetGene, MaxABC, CellTypesInCredibleSet,LeadVariant, !!sym(paste0("MaxABC.",cellgroup,".only")), !!sym(paste0(cellgroup,"InCredibleSet")))) %>% rename(CellTypesWithEPInteractions=CellTypesInCredibleSet) %>% rename(!!paste0(cellgroup,"CellTypesWithEPInteractions"):= !!paste0(cellgroup,"InCredibleSet"))
		ABC.df.core <- ABC.df.core %>% distinct() %>% mutate(TargetGene=ifelse(TargetGene=="MESDC1","TLNRD1",TargetGene))
		# Combine with ABC scores
		if (counter != 1){
	 		ABC.df.core <- ABC.df.core %>% select(-c(CellTypesWithEPInteractions, MaxABC))
		}
		CS.base.df <- left_join(CS.base.df, ABC.df.core, by=c("gene" = "TargetGene", "LeadVariant" = "LeadVariant")) 
		CS.base.df <- CS.base.df %>% select(-c(PeakChr,PeakStart,PeakEnd))
		CS.base.df <- CS.base.df %>% mutate(!!sym(paste0(cellgroup, "CellTypesWithEPInteractions")):=ifelse(!!sym(paste0(cellgroup,"CellTypesWithEPInteractions"))=="",NA,!!sym(paste0(cellgroup,"CellTypesWithEPInteractions")))) %>% mutate(!!sym(paste0("MaxABC.",cellgroup,".only")):=ifelse(!!sym(paste0("MaxABC.",cellgroup,".only"))==0,NA,!!sym(paste0("MaxABC.",cellgroup,".only"))))
	}
}

#############################################################
final.df <- CS.base.df
final.df  <- final.df %>% mutate(CredibleSet=paste0(Trait,"-",LeadVariant)) %>% relocate(CredibleSet, .before=LeadVariant) 
#############################################################
#############################################################
#Get the max abc score of each credible set
if ("MaxABC" %in% colnames(final.df)){
	final.df <- final.df %>% distinct() %>% group_by(LeadVariant,gene,Source) %>% mutate(MaxABC=ifelse(all(is.na(MaxABC)),NA,max(MaxABC, na.rm=TRUE))) %>% ungroup()
#Remove duplicates in CellTypes with EC Interactions
	final.df <- final.df %>% distinct() %>%  group_by(LeadVariant,gene,Source) %>% mutate(CellTypesWithEPInteractions = as.character(paste0(unique(CellTypesWithEPInteractions),collapse="|"))) %>% mutate(CellTypesWithEPInteractions = as.character(paste0(unique(unlist(strsplit(CellTypesWithEPInteractions, split="|",fixed = TRUE))),collapse = "|"))) %>% ungroup()
}
#For each cell group: 
for (cellgroup in cellgroups){
	if ((paste0("MaxABC.",cellgroup,".only") %in% colnames(final.df)) && (paste0(cellgroup,"CellTypesWithEPInteractions") %in% colnames(final.df)) && (paste0("MaxABC.",cellgroup,".only") %in% colnames(final.df))){
		#sometimes a gene maybe associated with multiple ABC scores because the same credible sets from different sources contain different number of variants. the following steps are to collapse these entries into one entry. 
		#this way it keeps the largest ABC score; The scores will be reranked at the end. 
		final.df <- final.df %>% distinct() %>% group_by(LeadVariant,gene,Source) %>% mutate(!!sym(paste0("MaxABC.",cellgroup,".only")):=ifelse(all(is.na(!!sym(paste0("MaxABC.",cellgroup,".only")))),NA,max(!!sym(paste0("MaxABC.",cellgroup,".only")), na.rm=TRUE))) %>% ungroup()
		#Again to collapse the cells with EP interactions. 
		final.df <- final.df %>% distinct() %>% group_by(LeadVariant,gene,Source) %>% mutate(!!sym(paste0(cellgroup,"CellTypesWithEPInteractions")):= as.character(paste0(unique(!!sym(paste0(cellgroup,"CellTypesWithEPInteractions"))),collapse="|"))) %>% mutate(!!sym(paste0(cellgroup,"CellTypesWithEPInteractions")) := as.character(paste0(unique(unlist(strsplit(!!sym(paste0(cellgroup,"CellTypesWithEPInteractions")), split="|",fixed = TRUE))),collapse = "|"))) %>% ungroup()
		#For each credible set, re-rank per gene max abc scores.
		final.df <- final.df %>% distinct() %>% group_by(LeadVariant,Source) %>% mutate(!!sym(paste0("MaxABC.Rank.",cellgroup,"Only")):=dense_rank(-!!sym(paste0("MaxABC.",cellgroup,".only")))) %>% ungroup() %>% relocate(!!sym(paste0("MaxABC.Rank.",cellgroup,"Only")), .after=!!sym(paste0("MaxABC.",cellgroup,".only")))
	}
}

#############################################################
# rank max abc
if ("MaxABC" %in% colnames(final.df)){
	final.df <- final.df %>% distinct() %>% group_by(LeadVariant,Source) %>% mutate(MaxABC.Rank=dense_rank(-MaxABC)) %>% ungroup() %>% relocate(MaxABC.Rank, .after=MaxABC)
}
# rank SNP to TSS
final.df <- final.df %>% arrange(chr, LeadVariantPos) %>% group_by(Source,LeadVariant) %>% mutate(rank_SNP_to_TSS=dense_rank(abs_SNP_to_TSS)) %>% ungroup() %>% as.data.frame()

#############################################################
write.table(final.df, file.path(opt$OutDir, paste0(opt$trait, "_cell2gene.txt")),sep="\t", quote=F, row.names=FALSE)

