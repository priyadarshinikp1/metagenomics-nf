#!/bin/bash
# generate_reports.sh
# Usage: ./generate_reports.sh downstream_results samples_down.csv
# Example: ./generate_reports.sh ./downstream_results samples_down.csv

set -euo pipefail

OUTDIR=$(realpath "$1")         # Base directory with all patient folders
SAMPLES=$(realpath "$2")        # CSV file listing patient IDs
RMD_TEMPLATE="/home/priyadarshini/metagenomics_nf/scripts/report_pdf.Rmd"

# Skip header, loop through each patient_id
tail -n +2 "$SAMPLES" | cut -d',' -f1 | while read -r PATIENT; do
    PATIENT_DIR="$OUTDIR/$PATIENT"

    # Check if patient folder exists
    if [ ! -d "$PATIENT_DIR" ]; then
        echo "⚠️  Missing analysis folder for $PATIENT, skipping"
        continue
    fi

    echo "📝 Generating report for $PATIENT ..."

    # Render R Markdown
    Rscript -e "rmarkdown::render(
        '$RMD_TEMPLATE',
        output_file='${PATIENT}_report.pdf',
        output_dir='$PATIENT_DIR',
        params=list(
            patient_name='$PATIENT',
            output_dir='$PATIENT_DIR'
        )
    )"

    echo "✅ Generated report for $PATIENT → $PATIENT_DIR/${PATIENT}_report.pdf"
done

