#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================
// Workflow: Database Setup
// =============================

process download_kneaddata {
    publishDir "${params.kneaddata_db}", mode: 'copy'

    output:
    path "*"

    script:
    """
    echo "[INFO] Downloading Kneaddata DB..."
    kneaddata_database --download ${params.kneaddata_db_type} bowtie2 .
    """
}

process download_kraken2 {
    publishDir "${params.kraken2_db}", mode: 'copy'

    output:
    path "*"

    script:
    """
    echo "[INFO] Downloading Kraken2 DB..."
    mkdir -p ${params.kraken2_db}
    wget -O kraken2_db.tar.gz ${params.kraken2_db_url}
    tar -xvzf kraken2_db.tar.gz -C ${params.kraken2_db}
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
    if [ -f "${params.kraken2_db}/database${params.bracken_read_len}mers.kmer_distrib" ]; then
        echo "[INFO] Bracken DB already exists. Skipping build."
        cp "${params.kraken2_db}/database${params.bracken_read_len}mers.kmer_distrib" .
    else
        echo "[INFO] Building Bracken DB..."
        bracken-build -d ${params.kraken2_db} -t ${task.cpus} -k ${params.bracken_kmer} -l ${params.bracken_read_len}
    fi
    """
}


// -----------------------------
// HUMAnN3 Databases
// -----------------------------
process download_humann_chocophlan {
    publishDir "${params.humann_chocophlan_db}", mode: 'copy'

    output:
    path "*"   // everything generated in work dir will be copied

    script:
    """
    echo "[INFO] Downloading HUMAnN3 ChocoPhlAn DB..."
    humann_databases --download chocophlan full .
    """
}


process download_humann_uniref {
    publishDir "${params.humann_uniref_db}", mode: 'copy'

    output:
    path "*"

    script:
    """
    echo "[INFO] Downloading HUMAnN3 UniRef90 DB..."
    humann_databases --download uniref ${params.humann_uniref_type} .
    """
}

// -----------------------------
// MetaPhlAn Database
// -----------------------------
process download_metaphlan_markers {
    
    script:
    """
    echo "[INFO] Downloading MetaPhlAn markers DB..."
    metaphlan --install .
    """
}

// =============================
// Workflow
// =============================
workflow {
    knead = download_kneaddata()
    kraken2_ch = download_kraken2()
    bracken = build_bracken(kraken2_ch)

    chocophlan = download_humann_chocophlan()
    uniref = download_humann_uniref()
    metaphlan = download_metaphlan_markers()
}
