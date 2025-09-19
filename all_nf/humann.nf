#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -------------------------------
// Parameters
// -------------------------------
params.samples_csv     = "./sample_humann.csv"
params.humann_out      = "./results/humann_out"
params.threads_humann  = 64

// -------------------------------
// Input samples CSV
// -------------------------------
samples_ch = Channel
    .fromPath(params.samples_csv)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, file(row.r1), file(row.r2)) }

// -------------------------------
// Merge paired reads
// -------------------------------
process merge_reads {
    cpus = 2
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq")

    script:
    """
    cat $r1 $r2 > ${sample_id}_merged.fastq
    """
}

// -------------------------------
// HUMAnN3 per-sample profiling
// -------------------------------
process run_humann {
    cpus = params.threads_humann
    tag "$sample_id"
    publishDir params.humann_out, mode: 'copy'

    input:
    tuple val(sample_id), path(merged)

    output:
    tuple val(sample_id),
          path("${sample_id}_merged_genefamilies.tsv"),
          path("${sample_id}_merged_pathabundance.tsv")

    script:
    """
    humann --input $merged --output ./ --threads ${task.cpus}
    """
}

// -------------------------------
// Post-processing (merge & normalize tables)
// -------------------------------
process humann_postprocessing {
    cpus = 4
    publishDir params.humann_out, mode: 'copy'

    input:
    path genefamilies_files
    path pathabundance_files
    val humann_out_dir 

    output:
    path "genefamilies.tsv"
    path "pathabundance.tsv"
    path "genefamilies_cpm.tsv"
    path "pathabundance_relab.tsv"

    script:
    """
    humann_join_tables --input $(dirname ${genefamilies_files[0]}) --output genefamilies.tsv --file_name genefamilies
    humann_join_tables --input $(dirname ${pathabundance_files[0]}) --output pathabundance.tsv --file_name pathabundance
    humann_renorm_table --input genefamilies.tsv --output genefamilies_cpm.tsv --units cpm
    humann_renorm_table --input pathabundance.tsv --output pathabundance_relab.tsv --units relab
   
    """
}

// -------------------------------
// Workflow
// -------------------------------

   workflow {
    merged_reads_ch = merge_reads(samples_ch)
    humann_results_ch = run_humann(merged_reads_ch)

    genefamilies_ch   = humann_results_ch.map { it[1] }
    pathabundance_ch  = humann_results_ch.map { it[2] }

    humann_postprocessing(genefamilies_ch, pathabundance_ch, Channel.value(params.humann_out))
}

