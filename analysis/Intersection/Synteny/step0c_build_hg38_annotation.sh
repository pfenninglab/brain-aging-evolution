#!/usr/bin/env bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hal

# ─────────────────────────────────────────────────────────────
# INPUTS
# ─────────────────────────────────────────────────────────────
DATADIR="/Aging/LQ_TACIT/data/"
#PEAKS_BED="${DATADIR}tables/Overlap_all_nuclei_neg_aging_pos.bed"
#PEAKS_BED="${DATADIR}tables/all_nuclei_bgd_427_random.bed"
#PEAKS_BED="${DATADIR}tables/all_nuclei_lq_neg_427_random.bed"
PEAKS_BED="${DATADIR}tables/all_nuclei_lq_neg_427_random2.bed"
#PEAKS_BED="/Aging/LQ_TACIT/phylolm/phylolm_r0_s1_bh_corrected_microglia_selected.bed"

TSS_BED="/evolution/VLTacit/data/raw_data/synteny_data/hg38_tss.bed"

OUTDIR="${DATADIR}liftover_out_random2/reference"
mkdir -p "$OUTDIR"

OUTFILE="${OUTDIR}/peakAnno_hg38.tsv"

# ─────────────────────────────────────────────────────────────
# 1. sort inputs (important for bedtools stability)
# ─────────────────────────────────────────────────────────────
PEAKS_SORTED="${OUTDIR}/hg38_peaks.sorted.bed"
TSS_SORTED="${OUTDIR}/hg38_tss.sorted.bed"

sort -k1,1 -k2,2n "$PEAKS_BED" > "$PEAKS_SORTED"
sort -k1,1 -k2,2n "$TSS_BED" > "$TSS_SORTED"

# ─────────────────────────────────────────────────────────────
# 2. nearest gene annotation
# ─────────────────────────────────────────────────────────────
bedtools closest \
    -a "$PEAKS_SORTED" \
    -b "$TSS_SORTED" \
    -D ref \
    -t first \
    > "${OUTDIR}/hg38.closest.bed"

# ─────────────────────────────────────────────────────────────
# 3. build peakAnno_df equivalent
# ─────────────────────────────────────────────────────────────
awk 'BEGIN{OFS="\t"}
{
    # assumes:
    # $4 = peak_id
    # $8 = gene symbol (from TSS annotation)
    # $9 = distance to TSS (from -D ref)
    print $4, "hg38", $8, $9
}' "${OUTDIR}/hg38.closest.bed" > "$OUTFILE"

# ─────────────────────────────────────────────────────────────
# 4. QC summary
# ─────────────────────────────────────────────────────────────
echo "DONE"
echo "Output: $OUTFILE"

echo "Total peaks:"
wc -l < "$PEAKS_BED"

echo "Annotated peaks:"
wc -l < "$OUTFILE"
