#!/bin/bash
# copy_all_patients.sh
# Usage: ./copy_all_patients.sh samples_down.csv

set -euo pipefail

SAMPLES=$1

# Base dir for FastQ Screen results
FASTQ_SCREEN_DIR="/home/microbiomeuser/stool_microbiome/results/fastq_screen"

# Skip header, loop through each patient_id
tail -n +2 "$SAMPLES" | cut -d',' -f1 | while read -r PATIENT; do
    OUTDIR="./downstream_results/$PATIENT"
    mkdir -p "$OUTDIR"

    # Source files
    KRAKEN_SRC="results/kraken2/${PATIENT}.kraken2.report"
    KRONA_PNG_SRC="results/krona/${PATIENT}_krona.png"
    KRONA_HTML_SRC="results/krona/${PATIENT}.krona.html"

    # FastQ Screen artifacts (PNG + HTML)
    FQS_P1_PNG="${FASTQ_SCREEN_DIR}/${PATIENT}_cleaned_paired_1_screen.png"
    FQS_P2_PNG="${FASTQ_SCREEN_DIR}/${PATIENT}_cleaned_paired_2_screen.png"
    FQS_P1_HTML="${FASTQ_SCREEN_DIR}/${PATIENT}_cleaned_paired_1_screen.html"
    FQS_P2_HTML="${FASTQ_SCREEN_DIR}/${PATIENT}_cleaned_paired_2_screen.html"

    # Copy Kraken2 report and convert to CSV with header
    if [[ -f "$KRAKEN_SRC" ]]; then
        cp "$KRAKEN_SRC" "$OUTDIR/"
        awk 'BEGIN {OFS=","; print "percentage,clade_count,direct_count,rank,taxid,name"} { $1=$1; print }' \
            "$KRAKEN_SRC" > "$OUTDIR/${PATIENT}_kraken2.csv"
    else
        echo "WARNING: Missing $KRAKEN_SRC"
    fi

    # Copy Krona artifacts
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

    # Copy FastQ Screen PNGs
    if [[ -f "$FQS_P1_PNG" ]]; then
        cp "$FQS_P1_PNG" "$OUTDIR/"
    else
        echo "WARNING: Missing $FQS_P1_PNG"
    fi

    if [[ -f "$FQS_P2_PNG" ]]; then
        cp "$FQS_P2_PNG" "$OUTDIR/"
    else
        echo "WARNING: Missing $FQS_P2_PNG"
    fi

    # Copy FastQ Screen HTMLs
    if [[ -f "$FQS_P1_HTML" ]]; then
        cp "$FQS_P1_HTML" "$OUTDIR/"
    else
        echo "WARNING: Missing $FQS_P1_HTML"
    fi

    if [[ -f "$FQS_P2_HTML" ]]; then
        cp "$FQS_P2_HTML" "$OUTDIR/"
    else
        echo "WARNING: Missing $FQS_P2_HTML"
    fi

    echo "? Processed $PATIENT ? $OUTDIR"
done
