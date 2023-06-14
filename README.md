# V2G
The V2G pipeline links genetic variants to their target genes on a cell-type specific basis. The only input required is a table of variants of interest, formatted as described in **Variant table** section.

You can download ABC predictions in 131 cell types and tissues from [here](https://www.engreitzlab.org/resources), and the corresponding accessible peaks from [here](https://mitra.stanford.edu/engreitz/public/SchnitzlerKang2023/EnhancerList.minus150). 

## Preprocessing 
Create a variant list for each trait.  This step is not required as long as the variant.list has the requreid columes specified in  **Variant table** section.  
*See **Variant table** section for required columns in the variant lists.*
LD-expand Aragam and Harst variants and include fine-mapped variants from the publication.
```
# LD-expand
bash preprocessing/CAD_*/log.addRsid.0.9.sh
# reformat
Rscript preprocessing/CAD_*/CreateVariantList.R
```
Preprocess the fine-mapped UK biobank variants.
```
# reformat
Rscript preprocessing/UKBB_*/CreateVariantList.R
```

## Run V2G
Create the environment from the snakmake/envs/V2G.yaml file
```
conda env create -f V2G.yaml
```
Run snakemake after setting up the config file
```
mkdir logs
snakemake --conda-frontend conda --profile sherlock --configfile snakemake/config/config.yaml --rerun-incomplete --snakefile snakemake/workflow/Snakefile --cluster "sbatch -n 1 -J {rule} -o logs/{rule}_{wildcards}.qout -e logs/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00"
```
### Config
```
VariantTable: describled in the **Variable table** section.
OutDir: the absolute path to the output file location.
CodeDir: the absolute path to the script directory. E.g. {path}/V2G/snakemake/workflow/scripts
hg19sizes: the absolute path to the chromosome sizes file. E.g. {path}/V2G/resources/hg19_chr_sizes
RemoveNonCoding: whether to remove non-coding genes. (TRUE/FALSE)

MungeCredibleSet:
        Genes: the absolute path to the bed file containing gene annotations. E.g., {path}/V2G/resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed
        Promoters: the absolute path to the file containing promoter annotations. E.g.,{path}/V2G/resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed"

UbiquitouslyExpressedGenes: the absolute path to the bed file containing ubiquitously expressed genes. E.g., {path}/V2G/resources/UbiquitouslyExpressedGenes.txt
ABCPRED: the absolute path to the ABC results for variant overlapping. E.g.{path}/V2G/resources/size_sorted_CombinedPredictions.AvgHiC.ABC0.015.minus150.txt.gz"

ALLPEAKS: the absolute path to the directory containing [cell type-specific peaks] (https://mitra.stanford.edu/engreitz/public/SchnitzlerKang2023/EnhancerList.minus150). 

CelltypeTable: the absolute path to the table specifying cell type name patterns for cell type groups. The table is for catagorizing the cell types in ABCPRD and ALLPEAKS. E.g. {path}/V2G/resources/grouped_celltype_table.txt"

LipidBloodAssociationTable: the absolute path to the table containing variants associated with lipid level regulation. E.g., {path}/V2G/resources/lipid.level.csv"

intersectWithABC:
  header: the absolute path to the file containing the header for the output file of the intersectWithABC step. E.g. {path}/V2G/resources/ABC_overlap.header"

intersectWithCelltypes:
        PIP: only consider variants with PIP larger than the value. E.g., 0

groupABCByCelltypes:
        PIP: same as above.

groupPeakOverlap:
        PIP: same as above
```

### Variant table
|ColumnName |Definition|Example|
|:---------:|:---------:|:------:|
|Trait      | The unique identifier of each trait. Required. |CAD_Aragam2021|
|fine_mapped_table | The absolute path to the fine-mapped variant list. Required. |variant.list.txt|
|LeadVariantCol | The lead variant column in the ***fine-mapped variant list***. Required. |LeadVariant|
|VariantCol | The variant ID column in the ***fine-mapped variant list***. Required. |RSID|
|ChrCol | The chromosome column in the ***fine-mapped variant list***. Required. |chr|
|PosCol | The variant position column in the ***fine-mapped variant list***. Required. |position|
|PCol |The p value column in the ***fine-mapped variant list***. Required. |P.value|
|PIPCol | The posterior probability column in the ***fine-mapped variant list***. NA if the column does not exist |P.value|
|Source | The source of the data |NA|
|ZeroIndexed | Whether the variant position is zero indexed (F/T). Required. |F|
|ExcludeVariants | IDs of variants to exclude from the analysis, if any. Required. |NA|


### R Session Info
```
R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /share/software/user/open/openblas/0.2.19/lib/libopenblasp-r0.2.19.so

Random number generation:
 RNG:     Mersenne-Twister
 Normal:  Inversion
 Sample:  Rejection

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] GenomicRanges_1.38.0 GenomeInfoDb_1.22.1  IRanges_2.20.2
 [4] S4Vectors_0.24.4     BiocGenerics_0.32.0  tibble_3.1.4
 [7] stringr_1.4.0        optparse_1.6.6       dplyr_1.0.7
[10] tidyr_1.1.2

loaded via a namespace (and not attached):
 [1] XVector_0.26.0         magrittr_2.0.1         zlibbioc_1.32.0
 [4] tidyselect_1.1.0       getopt_1.20.3          R6_2.5.1
 [7] rlang_0.4.11           fansi_0.5.0            tools_3.6.1
[10] utf8_1.2.2             DBI_1.1.0              ellipsis_0.3.2
[13] assertthat_0.2.1       lifecycle_1.0.1        crayon_1.4.1
[16] GenomeInfoDbData_1.2.2 purrr_0.3.4            bitops_1.0-7
[19] vctrs_0.3.8            RCurl_1.98-1.2         glue_1.6.2
[22] stringi_1.5.3          compiler_3.6.1         pillar_1.6.3
[25] generics_0.1.0         pkgconfig_2.0.3
```
### Main output files
- CAD_Aragam2021_cell2gene.txt: credible set to gene links. 
- Credibleset_gene_variant_info.tsv: variants in the credible set contributing to the credible set to gene links. 
- *relevant_variants.tsv: trait relevant variants in * cell type group.
- PeaksOverlapFull.Peaks.tsv: all variant-overlapping peaks.
- PeaksOverlapFull.tsv: all variants in peaks.
- CredibleSetPeakOverlapSummary.tsv: whether the credible sets overlapping with peaks in the cell type group. 
- CredibleSetPeakOverlapSummary.with*peakinfo.tsv: * cell type group peaks containing variants. 
- CredibleSetPeakOverlapSummary.withpeakinfo.only*.tsv: * cell type group peaks containing variants that only overlap peaks in this cell type group. 
- ABCOverlapFull.tsv: Vairant-overlapping ABC enhancers. 
- ABCOverlapFull.ranked.tsv: target genes of each credible set ranked by the maximum ABC scores of each gene in all cell types. 
- ABCOverlap.ranked.*.grouped.tsv: target genes of each credible set ranked by the maximum ABC scores of each gene in all cells in * cell type group. 
- ABCTopGeneTableWithCellCategories.tsv: cell type groups of top-ranked credible set-gene pairs (by ABC scores).
- ABCVariantOverlapSummary.tsv: a summary table of the ABC interactions of each variant.  
