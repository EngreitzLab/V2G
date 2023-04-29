suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
##############################################################################
project="Harst2017"
##############################################################################
## Load GWAS data and credible sets
variants <- read.delim("Harst2017_ld_expand.sumstats.tsv", header=F, stringsAsFactors=F)
colnames(variants) <- c("LeadVariant","RSID","chr","position","R2Leadvariant", "REF", "ALT", "Z", "P")

## There are some variants that are not found properly in Harst2021_ld_expand.sumstats.tsv (?)
variants <- subset(variants, !is.na(P))

variants <- variants %>%
  mutate(
    CredibleSet=paste0(project,"-",LeadVariant),
    LocusID=LeadVariant,
    Trait="CAD",
    variant=RSID,
    chr=paste0("chr",chr),
    position=as.numeric(position),
    P
    ) %>% distinct()

##############################################################################
write.table(variants, file="variant.list.txt", sep="\t", quote=F, row.names=F)
