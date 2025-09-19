#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samples   = "samples_down.csv"
params.outdir    = "downstream_results"
params.rscript   = "/home/priyadarshini/metagenomics_nf/scripts"
params.py        = "/home/priyadarshini/metagenomics_nf/scripts"
params.gutmgene  = "/home/priyadarshini/metagenomics_nf/databases/gutmgene_refs"
params.template  = "/home/priyadarshini/metagenomics_nf/scripts/report_pdf.Rmd"

/*************************************************************************
 * Input sample table -> channel
 *************************************************************************/
samples_ch = Channel
    .fromPath(params.samples)
    .splitCsv(header:true)
    .map { row ->
        tuple(
            row.patient_id,
            file(row.kraken2),
            file(row.bracken),
            file(row.genefamilies),
            file(row.control_kraken2),
            file(row.control_bracken),
            file(row.krona)
        )
    }

/*************************************************************************
 * 1) Run microbiome analysis
 *************************************************************************/
process microbiome_analysis {
    tag { patient_id }

    publishDir "${params.outdir}/${patient_id}", mode: 'copy'

    input:
    tuple val(patient_id), path(kraken2), path(bracken), path(genefamilies),
          path(control_kraken2), path(control_bracken), path(krona)

    output:
    tuple val(patient_id), path("${patient_id}"), emit: patient_dir

    script:
    """
    set -euo pipefail
    mkdir -p ${patient_id}
    Rscript ${params.rscript}/microbiome_analysis.R \
        --patient ${kraken2} \
        --bracken ${bracken} \
        --control ${control_kraken2} \
        --control_bracken ${control_bracken} \
        --genefamilies ${genefamilies} \
        --id ${patient_id} \
        --outdir ${patient_id}
    """
}

/*************************************************************************
 * 2) Rank network
 *************************************************************************/
process rank_network {
    tag { patient_id }

    publishDir "${params.outdir}/${patient_id}", mode: 'copy'

    input:
    tuple val(patient_id), path(rank_csv)

    output:
    path("${patient_id}_Microbe_Gene_rank_network.html")
    path("${patient_id}_Microbe_Gene_rank_network.png")

    script:
    """
    set -euo pipefail
    python ${params.py}/microbe_gene_rank_network.py
    """
}

/*************************************************************************
 * 3) MMGD mapping
 *************************************************************************/
process mmgd_mapping {
    tag { patient_id }

    publishDir "${params.outdir}/${patient_id}", mode: 'copy'

    input:
    tuple val(patient_id), path(relab_xlsx)

    output:
    tuple val(patient_id), path("${patient_id}_microbe_metabolite_gene_disease_network.csv"), emit: network
    tuple val(patient_id), path("${patient_id}_top_microbe_metabolite_gene_disease_network.html"), emit: html
    tuple val(patient_id), path("${patient_id}_top_microbe_metabolite_gene_disease_network.png"), emit: png

    script:
    """
    set -euo pipefail
    python ${params.py}/mmgd_fromR.py \
        ${relab_xlsx} \
        ${params.gutmgene}/Gut_Microbe_Microbial_metabolite.csv \
        ${params.gutmgene}/Microbial_metabolite_Host_Gene.csv \
        ${params.gutmgene}/Gut_Microbe_Host_Gene.csv
    """
}

/*************************************************************************
 * 4) MMGD downstream
 *************************************************************************/
process mmgd_downstream {

    tag { patient_id }

    publishDir "${params.outdir}/${patient_id}", mode: 'copy'

    input:
    tuple val(patient_id), path(network_csv)

    output:
    path("*"), emit: analysis  // Capture all outputs generated in this work dir

    script:
    """
    set -euo pipefail

    # Run downstream analysis
    python ${params.py}/mmgd_downstream.py

    # Move all generated outputs to current working dir
    # Assuming mmgd_downstream.py writes to 'analysis_outputs' subfolder
    if [ -d analysis_outputs ]; then
        mv analysis_outputs/* ./
        rmdir analysis_outputs
    fi
    """
}

/*************************************************************************
 * Workflow wiring
 *************************************************************************/
workflow {

    analysis_ch = microbiome_analysis(samples_ch)

    rank_net_input = analysis_ch.patient_dir.map { id, dir ->
        def f = file("${dir}/${id}_Microbe_Gene_RankSimilarity.csv")
        tuple(id, f)
    }.filter{ it[1].exists() }
    rank_network(rank_net_input)

    mmgd_input = analysis_ch.patient_dir.map { id, dir ->
        def f = file("${dir}/${id}_Microbiome_Filtered_RelativeAbundance.xlsx")
        tuple(id, f)
    }.filter{ it[1].exists() }
    mmgd_mapping(mmgd_input)

   // mmgd_downstream(mmgd_mapping.out.network)

}
