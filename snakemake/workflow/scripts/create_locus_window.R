suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--variants", type="character"),   
    make_option("--cs", type="character"),
    make_option("--outdir", type="character"),
    make_option("--trait", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))
############################################
all.cs.df <- read.table(opt$cs, sep = '\t', header=TRUE,stringsAsFactors=FALSE)
variants.df <- read.table(opt$variants, sep = '\t', header=TRUE,stringsAsFactors=FALSE) %>% select(CredibleSet, LeadVariant, LeadVariantPos) %>% distinct()

df <- all.cs.df %>% left_join(variants.df) %>% select(c(chr, CredibleSet, Trait,Source, LeadVariant, LeadVariantPos))

df <- df %>% mutate(region_start=LeadVariantPos-2000000) %>% mutate(region_end=LeadVariantPos+2000000) %>% relocate(region_start, .after=chr) %>% relocate(region_end,.after=region_start) %>% mutate(region_start=ifelse(region_start<0, 0, region_start))

write.table(df, file.path(opt$outdir, paste0(opt$trait, "_gwas_table_base.txt")),sep="\t", quote=F, row.names=FALSE)
df.bed <- df %>% select(chr, region_start, region_end, LeadVariant)
write.table(df.bed, file.path(opt$outdir, paste0(opt$trait, "_gwas_table_base.bed")), sep="\t", quote=F,row.names=FALSE, col.names=FALSE)



