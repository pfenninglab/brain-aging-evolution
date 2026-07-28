#!/bin/bash

# Get job id: e.g. "Submitted batch job 12345" -> "12345"
getjid() {
  echo "$1" | awk '{print $NF}'
}


echo "Phenotype is : $1"
PROJ_DIR=/projects/pfenninggroup/Aging/${1}_TACIT/
CODE_DIR=${PROJ_DIR}code/
OUT_DIR=${PROJ_DIR}phylolm/
PREDN_DIR=/projects/pfenninggroup/machineLearningForComputationalBiology/Cortex_Cell-TACIT/data/tidy_data/240_predictions_matrix_celltypes/
echo "$PREDN_DIR"

my_list=$(ls $PREDN_DIR)

for file in $my_list
do
    extracted_string=$(echo "$file")
    echo "$extracted_string"

    echo "1. Create sub-folder structure in the celltype output directory"
    mkdir ${OUT_DIR}${extracted_string}
    wait
    
    echo "2. Phylolm"
    cd ${CODE_DIR}

    echo "Round 1 -- Step 1 -- ocr_phylolm.r"
    jid1=$(sbatch -n 1 -p pool1 --wrap "source activate tacit; Rscript ${CODE_DIR}ocr_phylolm.r ${PROJ_DIR}data/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree ${PREDN_DIR}${extracted_string}/240_predictions_MatrixStacked.tsv ${PREDN_DIR}${extracted_string}/240_predictions_NamesList.txt ${PROJ_DIR}data/teeling_longevity_2-12-21.csv ${OUT_DIR}${extracted_string}/phylolm.csv 0 1 0 1 $1")
    echo "$jid1"

    echo "Round 1 -- Step 2 -- BH correction"
    jid1=$(getjid "$jid1")
    echo "$jid1"
    jid2=$(sbatch -n 1 -p pool1 --dependency=afterok:$jid1 --wrap "source activate tacit; Rscript bhCorrection_no_perm.R ${OUT_DIR}${extracted_string}/phylolm_r0_s1.csv ${OUT_DIR}${extracted_string}/phylolm_r0_s1_bh_corrected.csv")

    echo "Round 1 -- Step 3 — Prep for GREAT analysis"
    # Split positive and negative significant OCR after BH correction and convert to bed file
    jid2=$(getjid "$jid2")
    echo "$jid2"
    jid3=$(sbatch -n 1 -p pool1 --dependency=afterok:$jid2 --wrap "source activate tacit; python bh_corrected_ocr_bedfile_prep_no_perm.py ${1} ${extracted_string}")

done


