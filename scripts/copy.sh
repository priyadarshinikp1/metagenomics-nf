#!/bin/bash
# copy_all_patients.sh
# Usage: ./copy_all_patients.sh samples_down.csv

set -euo pipefail

SAMPLES=$1

# Skip header, loop through each patient_id
tail -n +2 "$SAMPLES" | cut -d',' -f1 | while read -r PATIENT; do
    OUTDIR="./downstream_results/$PATIENT"

    mkdir -p "$OUTDIR"

    KRAKEN_SRC="results/kraken2/${PATIENT}.kraken2.report"
    KRONA_PNG_SRC="results/krona/${PATIENT}_krona.png"
    KRONA_HTML_SRC="results/krona/${PATIENT}.krona.html"

    if [[ -f "$KRAKEN_SRC" ]]; then
        cp "$KRAKEN_SRC" "$OUTDIR/"
        awk 'BEGIN {OFS=","; print "percentage,clade_count,direct_count,rank,taxid,name"} { $1=$1; print }' \
            "$KRAKEN_SRC" > "$OUTDIR/${PATIENT}_kraken2.csv"
    else
        echo "WARNING: Missing $KRAKEN_SRC"
    fi

    if [[ -f "$KRONA_PNG_SRC" ]]; then
        cp "$KRONA_PNG_SRC" "$OUTDIR/"
    else
        echo "WARNING: Missing $KRONA_PNG_SRC"
    fi

    if [[ -f "$KRONA_HTML_SRC" ]]; then
        cp "$KRONA_HTML_SRC" "$OUTDIR/"
    else
        echo "WARNING: Missing $KRONA_HTML_SRC"
    fi

    echo "✅ Processed $PATIENT → $OUTDIR"
done

