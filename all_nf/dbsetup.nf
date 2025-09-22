#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================
// Workflow: Database Setup
// =============================

process download_kneaddata {
    publishDir "${params.kneaddata_db}", mode: 'copy'

    output:
    path "bowtie2*"

    script:
    """
    echo "[INFO] Checking Kneaddata DB..."
    if [ -d "${params.kneaddata_db}" ] && [ "\$(ls -A ${params.kneaddata_db})" ]; then
        echo "[INFO] Kneaddata DB already exists. Skipping download."
    else
        echo "[INFO] Downloading Kneaddata DB (${params.kneaddata_db_type})..."
        mkdir -p ${params.kneaddata_db}
        kneaddata_database --download ${params.kneaddata_db_type} bowtie2 ${params.kneaddata_db}
    fi
    """
}

process download_kraken2 {
    publishDir "${params.kraken2_db}", mode: 'copy'

    output:
    path "*"

    script:
    """
    echo "[INFO] Checking Kraken2 DB..."
    if [ -d "${params.kraken2_db}" ] && [ "\$(ls -A ${params.kraken2_db})" ]; then
        echo "[INFO] Kraken2 DB already exists. Skipping download."
    else
        echo "[INFO] Downloading Kraken2 DB..."
        mkdir -p ${params.kraken2_db}
        wget -O kraken2_db.tar.gz ${params.kraken2_db_url}
        tar -xvzf kraken2_db.tar.gz -C ${params.kraken2_db}
    fi
    """
}

process build_bracken {
    publishDir "${params.kraken2_db}", mode: 'copy'

    input:
    path kraken2_db

    output:
    path "database${params.bracken_read_len}mers.kmer_distrib"

    script:
    """
    echo "[INFO] Checking Bracken DB..."
    if [ -f "${params.kraken2_db}/database${params.bracken_read_len}mers.kmer_distrib" ]; then
        echo "[INFO] Bracken DB already exists. Skipping build."
    else
        echo "[INFO] Building Bracken DB..."
        bracken-build \
          -d ${params.kraken2_db} \
          -t ${task.cpus} \
          -k ${params.bracken_kmer} \
          -l ${params.bracken_read_len}
    fi
    """
}

// -----------------------------
// HUMAnN3 Databases
// -----------------------------
process download_humann_chocophlan {
    publishDir "${params.humann_chocophlan_db}", mode: 'copy'

    output:
    path "${params.humann_chocophlan_db}"

    script:
    """
    echo "[INFO] Checking HUMAnN3 ChocoPhlAn DB..."
    if [ -d "${params.humann_chocophlan_db}" ] && [ "\$(ls -A ${params.humann_chocophlan_db})" ]; then
        echo "[INFO] ChocoPhlAn DB already exists. Skipping download."
    else
        echo "[INFO] Downloading HUMAnN3 ChocoPhlAn DB..."
        mkdir -p ${params.humann_chocophlan_db}
        humann_databases --download chocophlan full ${params.humann_chocophlan_db}
    fi
    """
}

process download_humann_uniref {
    publishDir "${params.humann_uniref_db}", mode: 'copy'

    output:
    path "${params.humann_uniref_db}"

    script:
    """
    echo "[INFO] Checking HUMAnN3 UniRef90 DB..."
    if [ -d "${params.humann_uniref_db}" ] && [ "\$(ls -A ${params.humann_uniref_db})" ]; then
        echo "[INFO] UniRef90 DB already exists. Skipping download."
    else
        echo "[INFO] Downloading HUMAnN3 UniRef90 DB..."
        mkdir -p ${params.humann_uniref_db}
        humann_databases --download uniref ${params.humann_uniref_type} ${params.humann_uniref_db}
    fi
    """
}

// -----------------------------
// MetaPhlAn Database
// -----------------------------
process download_metaphlan_markers {
    publishDir "${params.metaphlan_db}", mode: 'copy'

    output:
    path "${params.metaphlan_db}"

    script:
    """
    echo "[INFO] Checking MetaPhlAn markers DB..."
    if [ -d "${params.metaphlan_db}" ] && [ "\$(ls -A ${params.metaphlan_db})" ]; then
        echo "[INFO] MetaPhlAn DB already exists. Skipping download."
    else
        echo "[INFO] Downloading MetaPhlAn markers DB..."
        mkdir -p ${params.metaphlan_db}
        metaphlan --install ${params.metaphlan_db}
    fi
    """
}

// =============================
// Workflow
// =============================
workflow {
    knead = download_kneaddata()
    kraken = download_kraken2()
    bracken = build_bracken(kraken)

    chocophlan = download_humann_chocophlan()
    uniref = download_humann_uniref()
    metaphlan = download_metaphlan_markers()
}
