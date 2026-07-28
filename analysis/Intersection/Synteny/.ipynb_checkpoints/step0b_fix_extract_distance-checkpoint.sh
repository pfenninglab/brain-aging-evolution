#!/usr/bin/env bash

DATADIR="/Aging/LQ_TACIT/data/"
OUTDIR="${DATADIR}liftover_out_random2"

for gdir in ${OUTDIR}/*/; do
genome=$(basename "$gdir")
closest="${gdir}closest.bed"
pct_file="${gdir}results.tsv"

#[[ -s "$pct_file" ]] || continue
pct=$(awk 'NR==1{print $NF}' "$pct_file")


[[ -s "$closest" ]] && \
awk -v g="$genome" -v pct="$pct" \
'BEGIN{OFS="\t"} {print g, $4, $8, $9, pct}' \
"$closest" > "${gdir}results.tsv"
done

cat ${OUTDIR}/*/results.tsv \
> ${OUTDIR}/closest_all.tsv

# Verify — should be 5 cols and a real distance in col 4
awk 'NR==1{print NF}' ${OUTDIR}/closest_all.tsv
awk 'NR<=3{print $1, $4}' ${OUTDIR}/closest_all.tsv
