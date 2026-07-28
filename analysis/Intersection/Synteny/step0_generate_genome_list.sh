#!/usr/bin/env bash
#SBATCH --partition=pfen3
#SBATCH --time=2-00:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%A.out

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hal

# ── paths ────────────────────────────────────────────────────────────────────
HAL="/projects/ikaplowlab/MetabolismTacit/cactus_alignments/vgp_577way_v1_alignment/vgp-577way-v1.hal"
REF_GENOME="GCA_000001405.15"

DATADIR="/Aging/LQ_TACIT/data/"
OUTDIR="${DATADIR}liftover_out_random2"
mkdir -p "$OUTDIR"

# ── 1. generate genome list once, then pick this array task's genome ──────────
GENOMES_FILE="${OUTDIR}/genomes.txt"
if [[ ! -f "$GENOMES_FILE" ]]; then
echo "Generating genome list..."

halStats --genomes "$HAL" | tr ' ' '\n' \
| grep -v "^$" | grep -v "^${REF_GENOME}$" | grep -v "Anc" \
> "$GENOMES_FILE"
fi



