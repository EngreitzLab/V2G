suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(GenomicRanges))
##############################################################################
project="Aragam2021"
##############################################################################
## Load GWAS data and credible sets
variants <- read.delim("Aragam2021_ld_expand.sumstats.tsv", header=F, stringsAsFactors=F)
header <- read.table(gzfile("FINAL_1MH_CAD_GWAS_summary_stats.tsv.gz"), header=F, nrows=1) %>% as.matrix() %>% as.character()
header <- c("LeadVariant","RSID","chr","position","R2Leadvariant", header)
colnames(variants) <- header

## There are some variants that are not found properly in Aragam2021_ld_expand.sumstats.tsv (?)
variants <- subset(variants, !is.na(N))

variants <- variants %>%
  mutate(
    CredibleSet=paste0(project,"-",LeadVariant),
    LocusID=LeadVariant,
    Disease="CAD",
    variant=RSID,
    chr=paste0("chr",chr),
    position=as.numeric(position),
    `P-value`
    ) %>% distinct() 

##############################################################################
write.table(variants, file="variant.list.txt",sep="\t", quote=F, row.names=F)
