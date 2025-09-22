# Shotgun Metagenomics Workflow (metagenomics-nf)

This repository contains a Nextflow-based pipeline for shotgun metagenomic analysis. It supports taxonomic profiling and functional profiling using HUMAnN3, with configuration files, sample manifests, environment settings and downstream steps.

---

## Contents

| Directory / File | Purpose |
|------------------|---------|
| `all_nf/` | Nextflow pipelines/scripts (e.g. `db_main.nf`, `dbsetup.nf`) |
| `data/` | Raw / intermediate data inputs (fastq, reference etc.) |
| `databases/` | Location for reference databases used by HUMAnN / Kraken etc. |
| `downstream_results/` | Processed output, summary tables, QC / visualization results etc. |
| `env/` | Environment files (e.g. `metagenomics_env.yml`) for reproducibility |
| `results/` | Main results from running pipelines (gene families, taxonomies etc.) |
| `scripts/` | Utility or helper scripts used in pipelines |
| `main.nf` | The main Nextflow workflow file controlling the pipeline |
| `down.nf` | Nextflow workflow for downstream analyses / plots |
| `nextflow.main.config` | Configuration for the main workflow (resources, inputs etc.) |
| `nextflow.down.config` | Config for downstream components |
| `samples.csv`, `sample_taxo.csv`, `sample_humann.csv`, `samples_down.csv` | Sample manifests / metadata files for different pipeline steps |

---

## Setup

1. **Clone this repo**

   ```bash
   git clone https://github.com/priyadarshinikp1/metagenomics-nf.git
   cd metagenomics-nf
   ```

2. **Install environment**

   There’s a conda environment under `env/` (presumably `metagenomics_env.yml`). Create and activate that:

   ```bash
   conda env create -f env/environment.yml 
   conda activate metagenomics_env
   ```

   Make sure all required tools (**Nextflow, HUMAnN3, Kraken2, Bracken, FastQC** etc.) are available via this environment.
   ```bash
   conda install nextflow
   ```

4. **Configure sample files and parameters**

   * Check / edit the sample manifest files (`samples.csv`, `sample_taxo.csv`, etc.) to reflect your samples, paths, metadata.
   * Update the config files (`nextflow.main.config`, `nextflow.down.config`) for your environment: directory paths, resources (threads, memory), database locations etc.
   * If needed, move or link your databases into the `databases/` directory (or point the configs to wherever they are stored).

---

## Running the Pipeline

Assuming you're using Nextflow:

### Main Workflow (Taxonomy + Functional Profiling)

```bash
nextflow run main.nf -c nextflow.main.config --samples samples.csv
```

You may need to pass additional params, for example manifests for taxonomic runs (`sample_taxo.csv`) or HUMAnN profiling (`sample_humann.csv`), depending on how `main.nf` is written.

### Downstream Analysis / Visualization

After the main workflow finishes:

```bash
nextflow run down.nf -c nextflow.down.config --samples samples_down.csv
```

This will likely generate summary tables, plots, reports, etc., based on outputs.

---

## Resume / Re-running

* Nextflow inherently supports resuming runs via the `-resume` flag:

  ```bash
  nextflow run main.nf -c nextflow.main.config --samples samples.csv -resume
  ```

* This will skip already completed tasks based on cached results.

---

## Output

You can expect your pipeline to produce (depending on which modules are enabled):

* **Taxonomic profiles** (e.g. classification outputs, abundance tables)
* **Functional profiles** (gene families, pathways via HUMAnN3)
* **Joined/normalized tables**
* **Quality control reports**
* **Visualizations** (plots, etc.)
* **Downstream analyses summaries**

These will reside under `results/` and/or `downstream_results/`, as per your config.

---

## Tips & Notes

* Verify sample manifest files carefully: correct sample IDs, file paths must match your data structure.
* Ensure you have enough computational resources (cores, memory) especially for heavy steps like functional profiling with large databases.
* Databases (HUMAnN, Kraken2 etc.) can be large; plan disk space.
* For reproducibility, record the versions of all software — Nextflow, tool versions in your conda env.

---

## Contact / Contribution

If you find bugs or want to contribute improvements (e.g. add new modules, adjust configs), feel free to open issues or pull requests.

---

## Example Workflow Diagram

```
Raw FASTQ ─► QC (FastQC, trimming) ─► Host filtering ─► Taxonomy (Kraken2/Bracken)  
               │                                    │  
               └─► Functional profiling (HUMAnN3) ◄──┘  
                     ▼  
            Joined / normalized tables + downstream summaries & plots  
```
