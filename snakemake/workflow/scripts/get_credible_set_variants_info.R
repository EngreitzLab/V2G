suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(optparse))
option.list <- list(
	make_option("--CS2Gene", type="character"),
	make_option("--variantList", type="character"),
    make_option("--peakOverlap", type="character"),
    make_option("--ABCOverlap", type="character"),
	make_option("--Cellgroups", type="character"),
	make_option("--OutDir", type="character"),
	make_option("--trait", type="character")
	)

opt <- parse_args(OptionParser(option_list=option.list))

cs2gene_df <- read.table(opt$CS2Gene, sep="\t", header=TRUE,stringsAsFactors=FALSE)
variant_table <- read.table(opt$variant, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cellgroups=as.list(strsplit(opt$Cellgroups, ",")[[1]])
peak_overlap_df <- read.table(opt$peakOverlap, sep="\t", header=TRUE, stringsAsFactors=FALSE) %>% select(variant,CredibleSet,CellType)
ABC_overlap_df <- read.table(opt$ABCOverlap, sep="\t", header=TRUE, stringsAsFactors=FALSE) %>% select(variant,CredibleSet,CellType, ABC.Score, TargetGene)
#############################################################
#filter variants for (coding|splice)
coding_df <- variant_table %>% filter(Coding) %>% select(variant,CredibleSet,Coding, CodingVariantGene) %>% mutate(codingVariant=ifelse(Coding, TRUE, FALSE)) %>% mutate_if(is.factor, as.character)

splice_df <- variant_table %>% filter(SpliceSite) %>% select(variant,CredibleSet, SpliceSite,SpliceSiteVariantGene) %>% mutate(splicingVariant=ifelse(SpliceSite, TRUE, FALSE)) %>% mutate_if(is.factor, as.character)
#############################################################
variant_peak_df <- cs2gene_df %>% select(CredibleSet,gene)
#Variant in Peaks
for (cellgroup in cellgroups){
    celltype_col=paste0("CellTypeOftheClosest", cellgroup, "Peak")
    variant_bool=paste0("variantInClosest", cellgroup, "Peak")
    df_peaks <- merge(cs2gene_df%>% select(CredibleSet, gene, !!sym(celltype_col)), peak_overlap_df, by.x=c("CredibleSet",celltype_col), by.y=c("CredibleSet","CellType")) %>% mutate(!!sym(variant_bool):=TRUE)
    variant_peak_df <- variant_peak_df %>% left_join(df_peaks) 
}

#############################################################
#Variant in ABC
df_ABC_interactions <- cs2gene_df %>% select(CredibleSet,gene, CellTypesWithEPInteractions) %>% separate_rows(CellTypesWithEPInteractions, sep="\\|", convert = TRUE) %>% as.data.frame() 
df_ABC_interactions <- merge(df_ABC_interactions, ABC_overlap_df, by.x=c("CredibleSet","CellTypesWithEPInteractions","gene"), by.y=c("CredibleSet", "CellType", "TargetGene")) %>% group_by(CredibleSet,gene) %>% mutate(ABC.rank = dense_rank(-ABC.Score)) %>% ungroup() %>% filter(ABC.rank==1) %>% select(-ABC.rank) %>% as.data.frame() %>% mutate(MaxABCVariant=TRUE)%>% group_by(CredibleSet) %>% mutate(MaxABC.Rank=dense_rank(-ABC.Score)) %>% rename(MaxABC.Score=ABC.Score) %>% ungroup() %>% as.data.frame() %>% mutate_if(is.factor, as.character)


variant_ABC_df <- cs2gene_df %>% select(CredibleSet,gene) 
for (cellgroup in cellgroups){
    celltype_col=paste0(cellgroup, "CellTypesWithEPInteractions")
    celltype_specific_df_ABC_interactions <- cs2gene_df %>% select(CredibleSet, gene, !!sym(celltype_col)) %>% separate_rows(., !!sym(celltype_col), sep="\\|", convert = TRUE) 
    celltype_specific_df_ABC_interactions <- merge(celltype_specific_df_ABC_interactions,ABC_overlap_df,by.x=c("CredibleSet",celltype_col, "gene"), by.y=c("CredibleSet", "CellType", "TargetGene")) %>% group_by(CredibleSet, gene) %>%  mutate(ABC.rank = dense_rank(-ABC.Score)) %>% ungroup() %>% filter(ABC.rank==1) %>% rename(!!paste0(cellgroup, ".ABC.Score"):=ABC.Score) %>% select(-ABC.rank) %>% as.data.frame()%>% mutate(!!sym(paste0("MacABC",cellgroup,"Variant")):=TRUE)
    variant_ABC_df <- variant_ABC_df %>% left_join(celltype_specific_df_ABC_interactions) %>% group_by(CredibleSet) %>% mutate(!!sym(paste0("MaxABC.Rank.", cellgroup, "Only")):=dense_rank(-!!sym(paste0(cellgroup, ".ABC.Score")))) %>% rename(!!paste0("MaxABC.Score.", cellgroup):=!!paste0(cellgroup, ".ABC.Score")) %>% ungroup() %>% as.data.frame() %>% mutate_if(is.factor, as.character)
}
#############################################################
if (nrow(coding_df)!=0 & nrow(splice_df)!=0){
    df <- cs2gene_df %>% select(CredibleSet,gene) %>% full_join(coding_df, by=c("CredibleSet", "gene"="CodingVariantGene")) %>% full_join(splice_df, by=c("CredibleSet", "gene"="SpliceSiteVariantGene", "variant")) %>% full_join(variant_peak_df,  by=c("CredibleSet","gene", "variant")) %>% full_join(df_ABC_interactions, by=c("CredibleSet", "gene", "variant")) %>% full_join(variant_ABC_df, by=c("CredibleSet", "gene", "variant")) %>% group_by(CredibleSet, gene, variant) %>% fill(everything(), .direction = "downup") %>% slice(1) %>% ungroup() %>% as.data.frame() %>% filter(!is.na(variant)) %>% ungroup() %>% as.data.frame()
} else if (nrow(coding_df)!=0 & nrow(splice_df)==0){
    df <- cs2gene_df %>% select(CredibleSet,gene) %>% full_join(coding_df, by=c("CredibleSet", "gene"="CodingVariantGene")) %>% full_join(variant_peak_df, by=c("CredibleSet","gene", "variant"))  %>% full_join(df_ABC_interactions, by=c("CredibleSet", "gene", "variant"))%>% full_join(variant_ABC_df, by=c("CredibleSet", "gene", "variant")) %>% group_by(CredibleSet, gene, variant) %>% fill(everything(), .direction = "downup") %>% slice(1) %>% ungroup() %>% as.data.frame() %>% filter(!is.na(variant)) %>% ungroup() %>% as.data.frame()
} else if  (nrow(coding_df)==0 & nrow(splice_df)!=0) {
    df <- cs2gene_df %>% select(CredibleSet,gene) %>% full_join(splice_df, by=c("CredibleSet", "gene"="SpliceSiteVariantGene")) %>% full_join(variant_peak_df, by=c("CredibleSet","gene", "variant")) %>% full_join(df_ABC_interactions, by=c("CredibleSet", "gene", "variant")) %>% full_join(variant_ABC_df, by=c("CredibleSet", "gene", "variant")) %>% group_by(CredibleSet, gene, variant) %>% fill(everything(), .direction = "downup") %>% slice(1) %>% ungroup() %>% as.data.frame() %>% filter(!is.na(variant)) %>% ungroup() %>% as.data.frame()
} else if (nrow(coding_df) !=0 & nrow(splice_df)!=0) {
    df <- cs2gene_df %>% select(CredibleSet,gene) %>% full_join(variant_peak_df,  by=c("CredibleSet","gene")) %>% full_join(df_ABC_interactions, by=c("CredibleSet", "gene", "variant"))  %>% full_join(variant_ABC_df, by=c("CredibleSet", "gene", "variant")) %>% group_by(CredibleSet, gene, variant) %>% fill(everything(), .direction = "downup") %>% slice(1) %>% ungroup() %>% as.data.frame() %>% filter(!is.na(variant)) %>% ungroup() %>% as.data.frame()
}

variant_table <- variant_table%>% select(chr, position, variant) %>%  distinct()

df_2 <- df %>% left_join(cs2gene_df %>% select(CredibleSet, LeadVariant, gene, SNP_to_TSS_distance,rank_SNP_to_TSS, starts_with("RankOfDistanceTo")),by=c("CredibleSet","gene")) %>% rename(lead_SNP_to_TSS_distance=SNP_to_TSS_distance) %>% rename(rank_lead_SNP_to_TSS=rank_SNP_to_TSS) %>% left_join(variant_table, by="variant") %>% relocate(position, .after=variant) %>% relocate(chr, .before=position) %>% distinct() %>% relocate(LeadVariant, .after=CredibleSet)

for (cellgroup in cellgroups){
    newcolName=paste0("RankOfGeneTSSDistanceTo", cellgroup, "PeakWithVariant")
    ogcolName=paste0("RankOfDistanceTo", cellgroup, "PeakWithVariant")
    df_2 <- df_2 %>% rename(!!newcolName :=!!ogcolName) %>% relocate(!!sym(newcolName), .after=!!sym(paste0("variantInClosest", cellgroup, "Peak")))
}
df_2 <- df_2 %>% distinct() %>% rename(TargetGene=gene)

write.table(df_2, file.path(opt$OutDir, "Credibleset_gene_variant_info.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

#############################################################
for (cellgroup in cellgroups){
    df_2_celltype <- df_2 %>% group_by(CredibleSet) %>% filter(!!sym(paste0("MaxABC.Rank.", cellgroup, "Only")) <= 2) %>% select(CredibleSet, LeadVariant, variant, position,!!sym(paste0(cellgroup, "CellTypesWithEPInteractions")), TargetGene, !!sym(paste0("MaxABC.Score.", cellgroup)),!!sym(paste0("MaxABC.Rank.", cellgroup, "Only")))
    write.table(df_2_celltype %>% distinct(), file.path(opt$OutDir, paste0(cellgroup, "_relevant_variants.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
}
