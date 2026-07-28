#!/usr/bin/env bash
#SBATCH --partition=pfen3
#SBATCH --time=2-00:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%A.out

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hal

DIR='/Aging/LQ_TACIT/data/'
DIR1='/Aging/LQ_TACIT/phylolm/'
DIR2='/Aging/Aging_SNAIL/data/threshold_98_2/'

# Aging Clock -- Differential OCR in vulnerable AD Inhibitory neurons (negative intrinsic age)
#bigWigAverageOverBed \
#  /Aging/LQ_TACIT/data/VGP/conservation/hg38/vgp-577way-v1-hg38.single.phyloP.ssrev.bw \
#  ${DIR2}positives/Anderson_PFC_AD_Inh_peaks_neg24_vs_all_hg38_positive.bed \
#  ${DIR}phylop_out/Anderson_PFC_AD_Inh_peaks_neg24_vs_all_hg38_positive_98_2_output.tsv

#bigWigAverageOverBed \
#  /Aging/LQ_TACIT/data/VGP/conservation/hg38/vgp-577way-v1-hg38.single.phyloP.ssrev.bw \
#  ${DIR2}negatives/Anderson_PFC_AD_Inh_peaks_neg24_vs_all_hg38_negative.bed \
#  ${DIR}phylop_out/Anderson_PFC_AD_Inh_peaks_neg24_vs_all_hg38_negative_98_2_output.tsv  

# LQ / Longevity -- Pooled associated OCR with LQ
#bigWigAverageOverBed \
#  /Aging/LQ_TACIT/data/VGP/conservation/hg38/vgp-577way-v1-hg38.single.phyloP.ssrev.bw \
#  ${DIR1}phylolm_r0_s1_bh_corrected_pos.bed \
#  ${DIR}phylop_out/phylolm_r0_s1_bh_corrected_pos_output.tsv

#bigWigAverageOverBed \
#  /Aging/LQ_TACIT/data/VGP/conservation/hg38/vgp-577way-v1-hg38.single.phyloP.ssrev.bw \
#  ${DIR1}phylolm_r0_s1_bh_corrected_neg.bed \
#  ${DIR}phylop_out/phylolm_r0_s1_bh_corrected_neg_output.tsv 

#-------------------------------------------------------------------------------------  

# Test peaks
#bigWigAverageOverBed \
  #/Aging/LQ_TACIT/data/VGP/conservation/hg38/vgp-577way-v1-hg38.single.phyloP.ssrev.bw \
  #${DIR}tables/Overlap_all_nuclei_neg_aging_pos.bed \
  #${DIR}phylop_out/Overlap_all_nuclei_neg_aging_pos_output.tsv
  
 
# Control peaks
bigWigAverageOverBed \
  /Aging/LQ_TACIT/data/VGP/conservation/hg38/vgp-577way-v1-hg38.single.phyloP.ssrev.bw \
  ${DIR}tables/all_nuclei_lq_neg_427_random2.bed \
  ${DIR}phylop_out/all_nuclei_lq_neg_427_random2_output.tsv # bgd 1
  #${DIR}phylop_out/all_nuclei_bgd_427_random_output.tsv # bgd 2
  #${DIR}phylop_out/longevity_pos.no_neg_overlap_output.tsv # bgd 3
  
  
  
