#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================
// Workflow: Database Setup
// =============================

process download_kneaddata {
    publishDir "${params.kneaddb}", mode: 'copy'

    output:
    path "bowtie2*"

    script:
    """
    echo "[INFO] Downloading Kneaddata DB (${params.kneaddb_type})..."
    mkdir -p ${params.kneaddb}
    kneaddata_database \
        --download ${params.kneaddb_type} bowtie2 ${params.kneaddb}
    """
}

process download_kraken2 {
    publishDir "${params.kraken2db}", mode: 'copy'

    output:
    path "*"

    script:
    """
    echo "[INFO] Downloading Kraken2 DB..."
    mkdir -p ${params.kraken2db}
    wget -O kraken2_db.tar.gz ${params.kraken2db_url}
    tar -xvzf kraken2_db.tar.gz -C ${params.kraken2db}
    """
}

process build_bracken {
    publishDir "${params.kraken2db}", mode: 'copy'

    input:
    path kraken2_db from download_kraken2.out

    output:
    path "database${params.bracken_read_len}mers.kmer_distrib"

    script:
    """
    echo "[INFO] Building Bracken DB..."
    bracken-build \
      -d ${params.kraken2db} \
      -t ${task.cpus} \
      -k ${params.bracken_kmer} \
      -l ${params.bracken_read_len}
    """
}

workflow {
    knead = download_kneaddata()
    kraken = download_kraken2()
    bracken = build_bracken(kraken)
}
