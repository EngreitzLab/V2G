suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))

option.list <- list(
	make_option("--outDir", type="character", help="Path to output directory"), 
	make_option("--grouped_celltype_table", type="character", help="Path to grouped cell type table"), 
	make_option("--cellGroup", type="character", help="CellCategory_to_group"), 
	make_option("--ranked_ABC_table", type="character", help="Path to ranked ABC table"), 
	make_option("--source", type="character", help="source if it's different from the one in credible set"), 
	make_option("--PIP", type="numeric", help="The PIP threshold to filter the variants"), 
	make_option("--helperFunctions", type="numeric")
)

opt <- parse_args(OptionParser(option_list=option.list))
setwd(opt$outDir)
source(opt$helperFunctions)
cellTags_df <- read.table(opt$grouped_celltype_table, header=TRUE, stringsAsFactors=F)
ABC.results <- read.table(opt$ranked_ABC_table, header=TRUE,stringsAsFactors=FALSE) 
if (!any(is.na(ABC.results$PosteriorProb))) {
	ABC.results <- ABC.results %>% filter(PosteriorProb >= opt$PIP)
}


ABC.grouped <- ABC.results %>% group_by(CredibleSet,TargetGene) %>% mutate(CellTypesInCredibleSet = as.character(paste0(sort(unique(CellType)), collapse = "|"))) %>% ungroup()

rownames(cellTags_df) <- cellTags_df$Cellgroup
celltype_tags <- unlist(strsplit(as.character(cellTags_df[opt$cellGroup, "Celltypes"]), split="|", fixed=TRUE))

ABC.grouped.cellgroup <- ABC.grouped %>% mutate(!!sym(opt$cellGroup):=ifelse(greplany(celltype_tags,CellType),1,0))

cell_group_max_ABC <- paste0("MaxABC.",opt$cellGroup, ".only")
celltypes_in_cell_group <- paste0(opt$cellGroup, "InCredibleSet")

# group by credible set, target gene, cell type group; calculate the largest ABC score and concatenate all the cell types
ABC.grouped.cellgroup <- ABC.grouped.cellgroup %>% group_by(CredibleSet,TargetGene,!!sym(opt$cellGroup)) %>% mutate(!!sym(cell_group_max_ABC):=max(ABC.Score)) %>% mutate(!!sym(celltypes_in_cell_group):=as.character(paste0(sort(unique(CellType)), collapse="|"))) %>% ungroup() 

# group by credible set and target gene,calculate the max cell type-specific ABC.
ABC.grouped.cellgroup <- ABC.grouped.cellgroup %>% group_by(CredibleSet, TargetGene) %>% mutate(!!sym(cell_group_max_ABC):=ifelse(!!sym(opt$cellGroup)==0, 0, !!sym(cell_group_max_ABC))) %>% mutate(!!sym(cell_group_max_ABC):=max(!!sym(cell_group_max_ABC))) %>% ungroup()

#Cell group specific cell types in the cell
ABC.grouped.cellgroup <-  ABC.grouped.cellgroup %>% group_by(CredibleSet, TargetGene) %>% mutate(!!sym(celltypes_in_cell_group):=ifelse(!!sym(opt$cellGroup) ==0, NA, !!sym(celltypes_in_cell_group))) %>% mutate(!!sym(celltypes_in_cell_group):=max(!!sym(celltypes_in_cell_group),na.rm = TRUE))%>% ungroup()

cell_group_max_ABC_rank <- paste0("MaxABC.Rank.",opt$cellGroup, ".only")
ABC.grouped.cellgroup <- ABC.grouped.cellgroup %>% group_by(CredibleSet) %>% mutate(!!sym(cell_group_max_ABC_rank):=dense_rank(-!!sym(cell_group_max_ABC))) %>% ungroup()

if (opt$source != ""){
 	ABC.grouped.cellgroup <- ABC.grouped.cellgroup %>% mutate(source=opt$source) %>% separate(CredibleSet, c("remove","LeadVariant"),sep="-", remove=FALSE) %>% select(-remove) %>% relocate(LeadVariant, .after=source)
 }else{
 	ABC.grouped.cellgroup <- ABC.grouped.celltype %>% separate(CredibleSet, c("source","LeadVariant"), sep="-", remove=FALSE)%>% relocate(LeadVariant, .after=source)
 }

output.table  <- paste0("ABCOverlap.ranked.", opt$cellGroup, ".grouped.tsv")
write.table(ABC.grouped.cellgroup, quote=FALSE, output.table, sep="\t", row.names=F)

