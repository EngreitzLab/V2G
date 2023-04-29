## Load GWAS data and credible sets

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
variants <- read.delim("DBP.bed", header=F)
variant.cols <- read.delim("UKBB_94traits_release1.cols", stringsAsFactors=F, header=F)$V1
colnames(variants) <- variant.cols

variants <- subset(variants, method == "SUSIE")
variants <- subset(variants, cs_id != -1)  ## These represent variants that SuSIE is not confident enough in to assign to a credible set

variants$LocusID <- variants$region
variants$CredibleSet <- paste0(variants$region,"-",variants$cs_id)
variants$Disease <- variants$trait
variants$variant <- variants$rsid
variants$chr <- variants$chromosome
variants$PosteriorProb <- variants$pip
variants$position <- variants$end
variants <- variants %>% group_by(CredibleSet) %>% mutate(rank=dense_rank(-PosteriorProb)) %>% mutate(LeadVariant=ifelse(rank==1, as.character(rsid), NA)) %>% mutate(LeadVariant=max(LeadVariant, na.rm=TRUE)) %>% ungroup() %>% select(-rank)

write.table(variants, file="variant.list.txt", quote=F, row.names=F, sep="\t")
