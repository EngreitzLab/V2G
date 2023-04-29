#!/bin/bash
REF=/oak/stanford/projects/genomics-refs/refs/
project="Aragam2021"
###############################
# Query 1000G for variants in LD with the lead variant:
LEADSNP=Aragam2021_Lead_SNP_info.txt
rm -f ${project}_ld_expand.tsv
rm -f ${project}_ld_expand_SNP_wo_ld.tsv
sed 1d ${LEADSNP} | \
cut -f 1,4 | sed 's/chr//' | \
while read CHR RSID; do
    echo ${RSID}
    plink --bfile ${REF}/ldscore/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR} \
          --r2 \
          --ld-snp ${RSID} \
          --ld-window-kb 1000 \
          --ld-window 99999 \
          --ld-window-r2 0.9 \
         --threads 1
    cat plink.ld | tr -s ' ' '\t' | sed 1d | sed 's/^\t//' | cut -f 3-7 >> ${project}_ld_expand.tsv
done

# Variants fine mapped by aragam; run the following to get the coordinates of the variants
ADDITIONAL_SNP="aragam_published_finemapped_variants.tsv"
sed 1d ${ADDITIONAL_SNP} | \
cut -f 1,3,4 | \
while read RSID CHR LeadSNP; do
	echo ${RSID}
	plink --bfile ${REF}/ldscore/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR} \
	      --r2 \
	      --ld-snp ${RSID} \
	      --ld-window-kb 1 \
	      --ld-window 99999 \
	      --ld-window-r2 1  \
	      --threads 1
    	if [ -f plink.ld ]; then
            echo -en ${LeadSNP} >> ${project}_ld_expand.tsv 
            #only keep the result that shows the coordiate of the variant itself
	    cat plink.ld | tr -s ' ' '\t' | sed 1d | sed 's/^\t//' | cut -f 3-7 | awk -F"\t" '$1 == $4 { print "\t"$2"\t"$3"\t"$4"\t"$5 }' | tr -d '\n' >> ${project}_ld_expand.tsv
	    echo -en "\n" >> ${project}_ld_expand.tsv
            rm plink.ld
    else
            echo $LeadSNP$'\t'""$'\t'""$'\t'""$'\t'"" >> ${project}_ld_expand_SNP_wo_ld.tsv
    fi
done
# Merge into summary stats table
cat ${project}_ld_expand.tsv | \
while read LEADSNP CHR POS SNP R2; do
    echo -en "${LEADSNP}\t${SNP}\t${CHR}\t${POS}\t${R2}\t"
    tabix FINAL_1MH_CAD_GWAS_summary_stats.tsv.gz ${CHR}:${POS}-${POS} | tr -d '\n'
    echo -en "\n"
done > ${project}_ld_expand.sumstats.tsv
