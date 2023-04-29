suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
option.list <- list(
  make_option("--helperFunctions", type="character"),
  make_option("--fineMappedVariants", type="character"),
  make_option("--genes", type="character"),
  make_option("--promoters", type="character"),
  make_option("--trait", type="character"),
  make_option("--outdir", type="character"),
  make_option("--leadVariantCol", type="character"),
  make_option("--variantCol", type="character"),
  make_option("--chrCol", type="character"),
  make_option("--positionCol", type="character"),
  make_option("--PCol", type="character"),
  make_option("--PIPCol", type="character"),
  make_option("--Source", type="character", help="Source of the data"),
  make_option("--zeroIndexed", type="logical", help="Whether the variant position is zero indexed. i.e. is the variant 1 bp off from the record in dbSNP"), 
  make_option("--excludeVariants", type="character", help="Variant to exclude")
)
opt <- parse_args(OptionParser(option_list=option.list))
source(opt$helperFunctions)
##############################################################################
variants=read.table(opt$fineMappedVariants, header=T, stringsAsFactors=F)
##############################################################################
## Load common data
genes <- readBed(opt$genes)
genes$symbol <- unlist(lapply(strsplit(as.character(as.matrix(genes$name)), ";"), "[", 1))
promoters <- readBed(opt$promoters)
promoters$symbol <- unlist(lapply(strsplit(as.character(as.matrix(promoters$name)), ";"), "[", 1))
promoters.gr <- GRangesFromBed(promoters)
############################################
names(variants)[names(variants) == opt$leadVariantCol] <- "LeadVariant"
names(variants)[names(variants) == opt$variantCol] <-  "variant"
names(variants)[names(variants) == opt$chrCol] <- "chr"
names(variants)[names(variants) == opt$positionCol] <- "position"
names(variants)[names(variants) == opt$PCol] <- "P"
variants=variants %>% subset(., select=which(!duplicated(names(.))))
if (opt$PIPCol != "None"){
	names(variants)[names(variants) == opt$PIPCol] <- "PosteriorProb"
} else {
	variants <- variants %>% mutate(PosteriorProb=NA)
}
variants=variants %>% subset(., select=which(!duplicated(names(.))))
print(head(variants))
if (opt$excludeVariants!="None"){
    excludeVariants=opt$excludeVariants %>% strsplit(split=",") %>% unlist()    
    variants <- variants %>% filter(!(LeadVariant %in% excludeVariants)) %>% filter(!is.na(P))
} else {
    variants <- variants %>% filter(!is.na(P))
}

if (opt$zeroIndexed){
    variants <- variants %>% mutate(
    CredibleSet=paste0(opt$trait, "-", LeadVariant),
    Trait=opt$trait,
    chr=ifelse(str_starts(chr, "Chr|chr"), chr, paste0("chr", chr)),
    position=as.numeric(position)+1) %>% group_by(LeadVariant) %>% mutate(LeadVariantPos_tmp=ifelse(variant==LeadVariant, position, NA)) %>% mutate(LeadVariantPos=max(as.numeric(LeadVariantPos_tmp), na.rm = TRUE)) %>% ungroup() %>% select(-LeadVariantPos_tmp)%>% filter(!is.infinite(LeadVariantPos))
} else {
    variants <- variants %>% mutate(
    CredibleSet=paste0(opt$trait,"-", LeadVariant),
    Trait=opt$trait,
    chr=ifelse(str_starts(chr, "Chr|chr"), chr, paste0("chr", chr)),
    position=as.numeric(position)) %>% group_by(LeadVariant) %>% mutate(LeadVariantPos_tmp=ifelse(variant==LeadVariant, position, NA)) %>% mutate(LeadVariantPos=max(as.numeric(LeadVariantPos_tmp), na.rm = TRUE)) %>% ungroup() %>% select(-LeadVariantPos_tmp) %>% filter(!is.infinite(LeadVariantPos))
}

variant.gr <- with(variants, GRangesFromBed(data.frame(chr=chr, start=position, end=position)))
variants2 <- annotateVariants(variants, variant.gr, promoters.gr, include.names=T)

variants <- variants2 %>% mutate(start=position, end=position, LocusID=LeadVariant) %>%
  dplyr:::select(
    chr,
    start, 
    end, 
    variant,
    P,
    PosteriorProb,
    CredibleSet,
    Trait,
    Coding,
    CodingVariantGene,
    SpliceSite,
    SpliceSiteVariantGene,
    Promoter,
    PromoterVariantGene,
    everything()
    ) 
##############################################################################
write.table(variants, file=file.path(opt$outdir,"variant.list.txt"), row.names=F, col.names=T, sep='\t', quote=F)

if (opt$PIPCol == "None")
variants <- variants2 %>% mutate(start=position, end=position, LocusID=LeadVariant) %>%
  dplyr:::select(
    chr,
    start,
    end,
    variant,
    P,
    PosteriorProb,
    CredibleSet,
    Trait,
    Coding,
    CodingVariantGene,
    SpliceSite,
    SpliceSiteVariantGene,
    Promoter,
    PromoterVariantGene,
    everything()
    ) %>% mutate(PosteriorProb=0)
all.cs <- makeCredibleSets(variants, genes, Source=opt$Source,include.names=TRUE)
write.table(all.cs, file=file.path(opt$outdir, "all.cs.txt"),  row.names=F, col.names=T, sep='\t', quote=F)
