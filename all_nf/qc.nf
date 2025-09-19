#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "./data/*_{R1,R2}.fastq.gz"
params.outdir = "./results"

// ----------------------
// Process 1: FastQC on raw reads
// ----------------------
process fastqc_raw {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "*.html", emit: html
    path "*.zip",  emit: zip

    script:
    """
    fastqc -t 4 -o ./ $read1 $read2
    """
}

// ----------------------
// Process 2: Fastp trimming
// ----------------------
process fastp {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.clean.fastq.gz"),
          path("${sample_id}_R2.clean.fastq.gz"),
          path("${sample_id}_fastp.json"),
          path("${sample_id}_fastp.html")

    script:
    """
    fastp \
        -i $read1 -I $read2 \
        -o ${sample_id}_R1.clean.fastq.gz \
        -O ${sample_id}_R2.clean.fastq.gz \
        --detect_adapter_for_pe \
        --cut_front --cut_tail --cut_mean_quality 20 \
        --length_required 50 \
        --trim_poly_g \
        --thread 4 \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    """
}

// ----------------------
// Process 3: FastQC on cleaned reads
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
    fastqc -t 4 -o ./ $r1 $r2
    """
}

// ----------------------
// Process 4: MultiQC
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
// Workflow
// ----------------------
workflow {

    // 1. Create paired-end channel from files
    reads_ch = Channel.fromFilePairs(params.reads, flat: true)

    // 2. FastQC on raw reads
    raw_reports = fastqc_raw(reads_ch)

    // 3. Fastp trimming
    cleaned = fastp(reads_ch)

    // 4. FastQC on cleaned reads
    clean_reports = fastqc_clean(cleaned.map { tuple(it[0], it[1], it[2]) }) // sample_id, R1, R2

    // 5. Combine FastQC zips + Fastp JSONs for MultiQC
    all_reports = raw_reports.zip.mix(clean_reports.zip).mix(cleaned.map { it[3] }) // index 3 is Fastp JSON

    // 6. Run MultiQC
    multiqc(all_reports)
}

