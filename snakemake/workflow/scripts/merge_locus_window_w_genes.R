suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--removeNonCoding", type="logical", default=TRUE, help="Whether to remove non-coding genes"),
    make_option("--helperFunctions", type="character"),
    make_option("--LocusWindow", type="character"),
    make_option("--GenesOverlappingLocusWindows", type="character"),
    make_option("--Outdir", type="character"), 
    make_option("--trait", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
source(opt$helperFunctions)

df <- read.table(opt$LocusWindow, sep = '\t', header=TRUE,stringsAsFactors=FALSE)
gene_overlap_table <- read.table(opt$GenesOverlappingLocusWindows, "\t", header=FALSE, stringsAsFactors=FALSE)
colnames(gene_overlap_table)=c("chr", "region_start", "region_end", "LeadVariant", "gene_chr", "pro_start", "pro_end","gene", "unknown", "strand")
if (opt$removeNonCoding){
	gene_overlap_table <- gene_overlap_table %>% filter(!greplany(c("^LINC","-AS","^MIR","RNU","^LOC"),gene)) %>% select(c(chr,region_start, region_end,LeadVariant, gene, pro_start, pro_end))
} else {
    gene_overlap_table <- gene_overlap_table %>% select(c(chr,region_start, region_end,LeadVariant, gene, pro_start, pro_end))
}
df_combined=left_join(df, gene_overlap_table,by=c("chr", "region_start", "region_end", "LeadVariant"))
df_combined <- df_combined %>% mutate(TSS_position=(pro_start+pro_end)/2) %>% mutate(SNP_to_TSS_distance=TSS_position-LeadVariantPos) %>% mutate(abs_SNP_to_TSS=abs(SNP_to_TSS_distance))
#keep at least two genes on each side of the lead variants (4 in total); include genes > +/- 500 kb away from the lead variants if there are < 2 genes on each side of the variant within the +/- 500kb window
df_combined <- df_combined %>% mutate(side=ifelse(SNP_to_TSS_distance<0, -1, 1)) %>% group_by(CredibleSet, side) %>% mutate(rank=dense_rank(abs_SNP_to_TSS)) %>% mutate(in_range=ifelse(abs_SNP_to_TSS<=500000, TRUE, FALSE)) %>% mutate(need_to_include_outside=ifelse((rank==1 | rank==2) & abs_SNP_to_TSS>500000, TRUE, FALSE)) %>% filter(in_range==TRUE|need_to_include_outside==TRUE) %>% mutate(side.sum=abs(sum(side))) %>% ungroup() %>% select(-c(side,rank,side.sum,in_range,need_to_include_outside))

df_combined <- df_combined %>% mutate(gene=ifelse(gene=="MESDC1","TLNRD1",gene)) %>% group_by(LeadVariant) %>% mutate(rank_SNP_to_TSS=dense_rank(abs_SNP_to_TSS))

write.table(df_combined, file.path(opt$Outdir, paste0(opt$trait, "_gwas_table_base_genes.txt")),sep="\t", quote=F, row.names=FALSE)
