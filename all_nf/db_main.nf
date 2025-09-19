#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------
// Parameters
// ----------------------
params.sample_csv        = "./samples.csv"
params.outdir            = "./results"
params.threads_qc        = 4
params.threads_host      = 4
params.threads_tax       = 64
params.kneaddata_db      = "/home/priyadarshinivizzy/metagenomics_nf/"
params.kneaddata_cleaned = "${params.outdir}/kneaddata_cleaned"
params.kraken2_db        = "/home/priyadarshinivizzy/metagenomics_nf/databases/kraken2"
params.bracken_read_len  = 150

// ----------------------
// Input Channel from CSV
// ----------------------
reads_ch = Channel
    .fromPath(params.sample_csv)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, file(row.fastq_r1), file(row.fastq_r2)) }

// ----------------------
// 1. FastQC on raw reads
// ----------------------
process fastqc_raw {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*.html", emit: html
    path "*.zip",  emit: zip

    script:
    """
    fastqc -t ${params.threads_qc} -o ./ $r1 $r2
    """
}

// ----------------------
// 2. Fastp trimming
// ----------------------
process fastp {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.clean.fastq.gz"),
          path("${sample_id}_R2.clean.fastq.gz"),
          path("${sample_id}_fastp.json"),
          path("${sample_id}_fastp.html")

    script:
    """
    fastp \
        -i $r1 -I $r2 \
        -o ${sample_id}_R1.clean.fastq.gz \
        -O ${sample_id}_R2.clean.fastq.gz \
        --detect_adapter_for_pe \
        --cut_front --cut_tail --cut_mean_quality 20 \
        --length_required 50 \
        --trim_poly_g \
        --thread ${params.threads_qc} \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    """
}

// ----------------------
// 3. FastQC on cleaned reads
// ----------------------
process fastqc_clean {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_clean", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*.html", emit: html
    path "*.zip",  emit: zip

    script:
    """
    fastqc -t ${params.threads_qc} -o ./ $r1 $r2
    """
}

// ----------------------
// 4. MultiQC
// ----------------------
process multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path reports

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc ${reports.join(' ')} -o ./ -n multiqc_report.html
    """
}

// ----------------------
// 5. Host Read Removal (Kneaddata)
// ----------------------
process kneaddata {
    tag "$sample_id"
    publishDir "${params.kneaddata_cleaned}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    val kneaddata_db //

    output:
    tuple val(sample_id),
          path("${sample_id}_cleaned_paired_1.fastq"),
          path("${sample_id}_cleaned_paired_2.fastq")

    script:
    """
    kneaddata \
        -i1 $r1 -i2 $r2 \
        -db $kneaddata_db \
        -o ./ \
        --output-prefix ${sample_id}_cleaned \
        -t ${params.threads_host} \
        -trimmomatic-threads 4

    """
}

// ----------------------
// 6. Kraken2
// ----------------------
process kraken2 {
    tag "$sample_id"
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    val kraken2_db

    output:
    tuple val(sample_id),
          path("${sample_id}.kraken2.out"),
          path("${sample_id}.kraken2.report")

    script:
    """
    kraken2 \
        --db $kraken2_db \
        --paired $r1 $r2 \
        --threads ${params.threads_tax} \
        --report ${sample_id}.kraken2.report \
        --output ${sample_id}.kraken2.out
    """
}

// ----------------------
// 7. Bracken
// ----------------------
process bracken {
    tag "$sample_id"
    publishDir "${params.outdir}/bracken", mode: 'copy'

    input:
    tuple val(sample_id), path(kraken_report)

    output:
    path "${sample_id}.bracken.species"

    script:
    """
    bracken \
        -d ${params.kraken2_db} \
        -i $kraken_report \
        -o ${sample_id}.bracken.species \
        -r ${params.bracken_read_len} -l S
    """
}

// ----------------------
// 8. Krona
// ----------------------
process krona {
    tag "$sample_id"
    publishDir "${params.outdir}/krona", mode: 'copy'

    input:
    tuple val(sample_id), path(kraken_out)

    output:
    path "${sample_id}.krona.html"

    script:
    """
    cut -f2,3 $kraken_out > ${sample_id}.krona.input
    ktImportTaxonomy ${sample_id}.krona.input -o ${sample_id}.krona.html
    """
}

// ----------------------
// Workflow
// ----------------------
workflow {

    // 1. QC: Raw FastQC
    raw_reports = fastqc_raw(reads_ch)

    // 2. Trimming
    trimmed = fastp(reads_ch)

    // 3. Clean FastQC
    clean_reports = fastqc_clean(trimmed.map { tuple(it[0], it[1], it[2]) })

    // 4. MultiQC: combine FastQC + Fastp JSON
    all_reports = Channel.empty()
    all_reports = all_reports.mix(raw_reports.html)
    all_reports = all_reports.mix(raw_reports.zip)
    all_reports = all_reports.mix(clean_reports.html)
    all_reports = all_reports.mix(clean_reports.zip)
    all_reports = all_reports.mix(trimmed.map { it[3] })  // Fastp JSON
    multiqc(all_reports)

    // 5. Host removal
    cleaned = kneaddata(trimmed.map { tuple(it[0], it[1], it[2]) }, params.kneaddata_db)

    // 6. Kraken2
    kraken_results = kraken2(cleaned, params.kraken2_db)

    // 7. Bracken
    bracken_input = kraken_results.map { sample_id, kraken_out, kraken_report -> tuple(sample_id, kraken_report) }
    bracken(bracken_input)

    // 8. Krona
    krona_input = kraken_results.map { sample_id, kraken_out, kraken_report -> tuple(sample_id, kraken_out) }
    krona(krona_input)
}

