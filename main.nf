#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------
// Input Channel from CSV
// ----------------------
reads_ch = Channel
    .fromPath(params.samples_csv)
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

    output:
    tuple val(sample_id),
          path("${sample_id}_cleaned_paired_1.fastq.gz"),
          path("${sample_id}_cleaned_paired_2.fastq.gz")

    script:
    """
    kneaddata \
        -i1 $r1 -i2 $r2 \
        -db ${params.kneaddata_db} \
        -o ./ \
        --output-prefix ${sample_id}_cleaned \
        -t ${params.threads_host} \
        --remove-intermediate-output \
        --verbose

    gzip ${sample_id}_cleaned_paired_1.fastq
    gzip ${sample_id}_cleaned_paired_2.fastq
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

    output:
    tuple val(sample_id),
          path("${sample_id}.kraken2.out"),
          path("${sample_id}.kraken2.report")

    script:
    """
    kraken2 \
        --db ${params.kraken2_db} \
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
// 9. Merge reads for HUMAnN3
// ----------------------
process merge_reads {
    cpus = 2
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq")

    script:
    """
    zcat $r1 $r2 > ${sample_id}_merged.fastq
    """
}

// ----------------------
// 10. HUMAnN3 profiling
// ----------------------
process run_humann {
    cpus = params.threads_humann
    tag "$sample_id"
    publishDir params.humann_out, mode: 'copy'

    input:
    tuple val(sample_id), path(merged)

    output:
    tuple path("${sample_id}_genefamilies.tsv"),
          path("${sample_id}_pathabundance.tsv")

    script:

    """
    humann --input $merged --output ./ --threads ${task.cpus} --metaphlan-options "--bowtie2db ${params.metaphlan_db}"
    mv ${sample_id}_merged_genefamilies.tsv ${sample_id}_genefamilies.tsv
    mv ${sample_id}_merged_pathabundance.tsv ${sample_id}_pathabundance.tsv
    rm $merged
    """
}

// ----------------------
// 11. HUMAnN3 postprocessing
// ----------------------
process humann_postprocessing {
    cpus = 4
    publishDir params.humann_out, mode: 'copy'

    input:
    path genefamilies_files

    output:
    path "genefamilies.tsv"
    path "pathabundance.tsv"
    path "genefamilies_cpm.tsv"
    path "pathabundance_relab.tsv"

    script:
    """
    humann_join_tables --input /home/priyadarshini/metagenomics_nf/res/humann_out --output genefamilies.tsv --file_name genefamilies
    humann_join_tables --input /home/priyadarshini/metagenomics_nf/res/humann_out --output pathabundance.tsv --file_name pathabundance
    humann_renorm_table --input genefamilies.tsv --output genefamilies_cpm.tsv --units cpm
    humann_renorm_table --input pathabundance.tsv --output pathabundance_relab.tsv --units relab
    """
}

// ----------------------
// Workflow
// ----------------------
workflow {
    // 1. QC
    raw_reports   = fastqc_raw(reads_ch)
    trimmed       = fastp(reads_ch)
    clean_reports = fastqc_clean(trimmed.map { tuple(it[0], it[1], it[2]) })
    all_reports   = raw_reports.zip
                     .mix(clean_reports.zip)
                     .mix(trimmed.map { it[3] })
    multiqc(all_reports)

    // 2. Host removal
    cleaned = kneaddata(trimmed.map { tuple(it[0], it[1], it[2]) })

    // 3. Taxonomy profiling
    kraken_results = kraken2(cleaned)
    bracken_input  = kraken_results.map { sample_id, kraken_out, kraken_report -> tuple(sample_id, kraken_report) }
    bracken(bracken_input)
    krona_input    = kraken_results.map { sample_id, kraken_out, kraken_report -> tuple(sample_id, kraken_out) }
    krona(krona_input)

    // 4. Functional profiling
    merged_reads   = merge_reads(cleaned)
    humann_results = run_humann(merged_reads)

    // Postprocess HUMAnN3 results
    humann_postprocessing(humann_results.collect())
}

