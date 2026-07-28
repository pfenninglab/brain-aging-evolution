#!/bin/bash
source activate tacit

DATA_DIR=/Aging/Data/Anderson/Multiomics/ATAC/ArchR/Output/
#my_list=$(ls "$DATA_DIR"/PeakSet*.bed 2>/dev/null)
my_list=$(ls "$DATA_DIR"/Anderson_PeakMatrix_Reads.bed 2>/dev/null)

anno=/Aging/Data/Roadmap_Epigenomics/ChromHMM_E073_BRN.DL.PRFRNTL.CRTX_hg38lift/E073_18_core_K27ac_hg38lift_mnemonics.bed.gz

for file in $my_list
do
    # Get filename without path and extension
    filename=$(basename "$file" .bed)
    echo "$filename"
    
    bedtools intersect -a "$file" -b "$anno" -wb | \
        awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7}' > "${DATA_DIR}/${filename}_rm_epigenomics_anno_intersect.txt"
        
    bedtools intersect -a "$file" -b "$anno" -c  > "${DATA_DIR}/${filename}_rm_epigenomics_anno_intersect_counts.txt"
        
done

