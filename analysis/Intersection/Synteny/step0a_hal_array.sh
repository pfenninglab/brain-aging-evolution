#!/usr/bin/env bash
#SBATCH --partition=pfen3
#SBATCH --time=2-00:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-576%30
#SBATCH --output=logs/%A_%a.out

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hal

# ── paths ────────────────────────────────────────────────────────────────────
HAL="/MetabolismTacit/cactus_alignments/vgp_577way_v1_alignment/vgp-577way-v1.hal"
REF_GENOME="GCA_000001405.15"

DATADIR="/Aging/LQ_TACIT/data/"
#PEAKS_BED="${DATADIR}tables/Overlap_all_nuclei_neg_aging_pos.bed"
#PEAKS_BED="${DATADIR}tables/all_nuclei_bgd_427_random.bed"
#PEAKS_BED="${DATADIR}tables/all_nuclei_lq_neg_427_random.bed"
PEAKS_BED="${DATADIR}tables/all_nuclei_lq_neg_427_random2.bed"
#PEAKS_BED="/Aging/LQ_TACIT/phylolm/phylolm_r0_s1_bh_corrected_microglia_selected.bed"

TSS_BED="/evolution/VLTacit/data/raw_data/synteny_data/hg38_tss.bed"

OUTDIR="${DATADIR}liftover_out_random2"
#mkdir -p "$OUTDIR"

# ── 1. generate genome list once, then pick this array task's genome ──────────
GENOMES_FILE="${OUTDIR}/genomes.txt"

GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GENOMES_FILE")
echo "Array task ${SLURM_ARRAY_TASK_ID}: processing ${GENOME}"

N_PEAKS=$(wc -l < "$PEAKS_BED")

# ── 2. per-genome work ────────────────────────────────────────────────────────
GDIR="${OUTDIR}/${GENOME}"
mkdir -p "$GDIR"
GENOME_RESULTS="${GDIR}/results.tsv"

# 2a. lift enhancer peaks
PEAKS_LIFTED="${GDIR}/peaks_lifted.bed"
halLiftover "$HAL" "$REF_GENOME" "$PEAKS_BED" "$GENOME" "$PEAKS_LIFTED" \
2>"${GDIR}/peaks_liftover.log" || true

# 2b. lift TSS annotations
TSS_LIFTED="${GDIR}/tss_lifted.bed"
halLiftover "$HAL" "$REF_GENOME" "$TSS_BED" "$GENOME" "$TSS_LIFTED" \
2>"${GDIR}/tss_liftover.log" || true

# compute lift rate (even if 0)
N_LIFTED=0
[[ -s "$PEAKS_LIFTED" ]] && N_LIFTED=$(wc -l < "$PEAKS_LIFTED")
PCTLIFTED=$(echo "scale=2; $N_LIFTED / $N_PEAKS * 100" | bc)
echo "  ${N_LIFTED}/${N_PEAKS} peaks lifted (${PCTLIFTED}%)"

# if either liftover produced nothing, write sentinel and exit
if [[ ! -s "$PEAKS_LIFTED" || ! -s "$TSS_LIFTED" ]]; then
echo "  NOTE: skipping closest-gene step for ${GENOME} (no lifted intervals)"
echo -e "${GENOME}\tNA\tNA\tNA\t${PCTLIFTED}" > "$GENOME_RESULTS"
exit 0
fi

# 2c. sort both files for bedtools
PEAKS_SORTED="${GDIR}/peaks_sorted.bed"
TSS_SORTED="${GDIR}/tss_sorted.bed"
sort -k1,1 -k2,2n "$PEAKS_LIFTED" > "$PEAKS_SORTED"
sort -k1,1 -k2,2n "$TSS_LIFTED"   > "$TSS_SORTED"

# 2d. find nearest gene to each lifted peak
CLOSEST="${GDIR}/closest.bed"
bedtools closest \
-a "$PEAKS_SORTED" \
-b "$TSS_SORTED" \
-D ref \
-t first \
> "$CLOSEST"

# 2e. write per-genome results file
awk -v genome="$GENOME" -v pct="$PCTLIFTED" \
'BEGIN{OFS="\t"} {print genome, $4, $8, pct}' \
"$CLOSEST" > "$GENOME_RESULTS"

echo "  Done: results written to ${GENOME_RESULTS}"

#---------------------------------------------------------------

