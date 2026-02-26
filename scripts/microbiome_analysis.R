
# =========================================================
# 1. Load Required Libraries
# =========================================================
library(tidyverse)
library(data.table)
library(writexl)
library(readxl)
library(igraph)
library(visNetwork)
library(htmlwidgets)
library(readr)
library(stringr)
library(vegan)
library(pheatmap)
library(rmarkdown)
library(knitr)
library(patchwork)
library(kableExtra)

# =========================================================
# 2. Parse Command-Line Arguments
# =========================================================
args <- commandArgs(trailingOnly = TRUE)

patient_name <- NULL
patient_file <- NULL
control_file <- NULL
genefamilies_file <- NULL
outdir <- NULL
bracken_file <- NULL
control_bracken_file <- NULL

for (i in seq(1, length(args), by=2)) {
  if (args[i] == "--id") patient_name <- args[i+1]
  else if (args[i] == "--patient") patient_file <- args[i+1]
  else if (args[i] == "--control") control_file <- args[i+1]
  else if (args[i] == "--genefamilies") genefamilies_file <- args[i+1]
  else if (args[i] == "--outdir") outdir <- args[i+1]
  else if (args[i] == "--bracken") bracken_file <- args[i+1]
  else if (args[i] == "--control_bracken") control_bracken_file <- args[i+1]
}

if (is.null(patient_name) || is.null(patient_file) || is.null(control_file) || is.null(genefamilies_file)) {
  stop("Required arguments missing: --id, --patient, --control, --genefamilies")
}
if (is.null(outdir)) outdir <- "results"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# =========================================================
# 3. Function to Load Kraken2 Reports
# =========================================================
load_kraken_report <- function(filepath, sample_name) {
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))
  df <- fread(filepath, header = FALSE, sep = "\t")
  colnames(df) <- c("Percent", "Reads_Clade", "Reads_Direct", "Rank", "TaxID", "Name")
  df$Sample <- sample_name
  return(df)
}

sample_df <- load_kraken_report(patient_file, "Patient")
control_df <- load_kraken_report(control_file, "Control")

# =========================================================
# 4. Combine, Filter for Species-Level & Compute log2FC
# =========================================================
species_df <- bind_rows(sample_df, control_df) %>%
  #filter(Rank == "S") %>%
  filter(Rank %in% c("S", "G")) %>%
  select(Sample, Name, Reads_Clade) %>%
  pivot_wider(names_from = Sample, values_from = Reads_Clade, values_fill = 0)

species_counts <- species_df %>%
  column_to_rownames("Name") %>%
  mutate(log2FoldChange = log2((Patient + 1) / (Control + 1))) %>%
  arrange(desc(abs(log2FoldChange)))

significant_microbes <- species_counts %>%
  filter(log2FoldChange > 2) %>%
  rownames_to_column("Name")

patient_diff_abundant <- significant_microbes %>%
  left_join(sample_df, by="Name") %>%
  select(Name, Patient, Control, log2FoldChange, Percent, Reads_Clade, Reads_Direct, Rank, TaxID)

write_xlsx(patient_diff_abundant, file.path(outdir, paste0(patient_name, "_Differentially_Abundant_Microbes.xlsx")))


# =========================================================
# Heatmap (Patient vs Control)
# Use log10(count+1); show top 30 by total abundance.
# =========================================================
top_hm <- 30
hm_df <- species_df %>%
  mutate(Total = Patient + Control) %>%
  arrange(desc(Total)) %>%
  slice_head(n = top_hm) %>%
  select(Name, Patient, Control) %>%
  column_to_rownames("Name")

# Transform for visualization
hm_mat <- log10(as.matrix(hm_df) + 1)

pheatmap(hm_mat,
         filename = file.path(outdir, paste0(patient_name, "_heatmap_top", top_hm, ".png")),
         width = 6, height = 8, dpi = 300,
         main = paste0("Top ", top_hm, " taxa (log10 counts)"),
         cluster_rows = TRUE, cluster_cols = FALSE,
         color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100))
# =========================================================
# 5. Top Species Analysis
# =========================================================
species_top <- species_df %>%
  group_by(Name) %>%
  summarise(PatientReads=sum(Patient), ControlReads=sum(Control), .groups="drop") %>%
  arrange(desc(PatientReads)) %>%
  slice_head(n=10)

# Top 10 patient species
p1_file <- file.path(outdir, paste0(patient_name, "_Top10_Patient_Species.png"))
p1 <- ggplot(species_top, aes(x=reorder(Name, PatientReads), y=PatientReads)) +
  geom_bar(stat="identity", fill="steelblue") + coord_flip() +
  labs(title="Top 10 Species in Patient Sample", x="Species", y="Reads") + theme_minimal()
ggsave(p1_file, p1, width=8, height=6, dpi=300)

# Top 10 patient vs control
species_long <- species_top %>%
  pivot_longer(cols=c(PatientReads, ControlReads), names_to="Sample", values_to="Reads")
p2_file <- file.path(outdir, paste0(patient_name, "_Top10_Patient_vs_Control.png"))
p2 <- ggplot(species_long, aes(x=reorder(Name, Reads), y=Reads, fill=Sample)) +
  geom_bar(stat="identity", position="dodge") + coord_flip() +
  labs(title="Top 10 Species: Patient vs Control", x="Species", y="Reads") + theme_minimal()
ggsave(p2_file, p2, width=8, height=6, dpi=300)

# =========================================================
# 6. Relative Abundance
# =========================================================
microbiome_filtered <- patient_diff_abundant %>%
  filter(Percent>0) %>%
  mutate(RelativeAbundance = Reads_Clade / sum(Reads_Clade) * 100)
write_xlsx(microbiome_filtered, file.path(outdir, paste0(patient_name, "_Microbiome_Filtered_RelativeAbundance.xlsx")))

# =========================================================
# 7. Load HUMAnN Gene Families
# =========================================================
gf_cleaned <- read_tsv(genefamilies_file, comment="#", col_names=c("Gene_Family","Abundance"))

gf_strat <- gf_cleaned %>%
  filter(str_detect(Gene_Family,"\\|")) %>%
  separate(Gene_Family, into=c("GeneID","Microbe"), sep="\\|") %>%
  rename(GeneAbundance=Abundance) %>%
  mutate(GeneID=str_replace(GeneID,"^UniRef\\d+_",""),
         Microbe_clean=Microbe %>% str_replace("g__","") %>% str_replace("s__","") %>%
           str_replace_all("_"," ") %>% str_replace("^[^.]+\\.",""))

# =========================================================
# 8. Rank-based Integration
# =========================================================
microbe_ranks <- microbiome_filtered %>% mutate(MicrobeRank=rank(-RelativeAbundance)) %>% select(Name, MicrobeRank)
gene_ranks <- gf_strat %>% group_by(Microbe_clean) %>% mutate(GeneRank=rank(-GeneAbundance)) %>% ungroup()
comparison_rank <- gene_ranks %>% inner_join(microbe_ranks, by=c("Microbe_clean"="Name"))
rank_similarity <- comparison_rank %>% mutate(RankDifference=abs(MicrobeRank-GeneRank),
                                              RankSimilarity=1-(2*RankDifference/max(RankDifference,na.rm=TRUE))) %>%
  select(Microbe_clean,GeneID,MicrobeRank,GeneRank,RankDifference,RankSimilarity) %>% arrange(RankDifference)
write_csv(rank_similarity, file.path(outdir, paste0(patient_name, "_Microbe_Gene_RankSimilarity.csv")))

# =========================================================
# 9. Beta Diversity (Bray-Curtis & Jaccard)
# =========================================================
if (is.null(bracken_file)) {
  bracken_file <- str_replace(patient_file,".kraken2.report$",".bracken.species")
}
if (is.null(control_bracken_file)) {
  control_bracken_file <- str_replace(control_file,".kraken2.report$",".bracken.species")
}

patient_bracken <- read_tsv(bracken_file)
control_bracken <- read_tsv(control_bracken_file)

patient_df <- patient_bracken %>% select(name,new_est_reads) %>% rename(Patient=new_est_reads)
control_df <- control_bracken %>% select(name,new_est_reads) %>% rename(Control=new_est_reads)
species_matrix <- full_join(patient_df, control_df, by="name") %>% replace_na(list(Patient=0,Control=0)) %>% column_to_rownames("name") %>% t()

bray <- vegdist(species_matrix, method="bray")
jacc <- vegdist(species_matrix, method="jaccard", binary=TRUE)

write_xlsx(list(
  BrayCurtis=as.data.frame(as.matrix(bray)) %>% rownames_to_column("Sample"),
  Jaccard=as.data.frame(as.matrix(jacc)) %>% rownames_to_column("Sample")
), file.path(outdir, paste0(patient_name,"_BetaDiversity.xlsx")))

# =========================================================
# Beta Diversity (Bray-Curtis) using absolute differences
# =========================================================

# Prepare species matrix for beta diversity
species_df_beta <- species_df %>%
  select(Name, PatientReads = Patient, ControlReads = Control) %>%
  column_to_rownames("Name")

# Transpose for vegan
species_matrix <- t(species_df_beta)

# Compute Bray-Curtis distance
bray <- vegdist(species_matrix, method = "bray")
jacc <- vegdist(species_matrix, method = "jaccard", binary = TRUE)

# Save distance matrices to Excel
write_xlsx(list(
  BrayCurtis = as.data.frame(as.matrix(bray)) %>% rownames_to_column("Sample"),
  Jaccard = as.data.frame(as.matrix(jacc)) %>% rownames_to_column("Sample")
), file.path(outdir, paste0(patient_name,"_BetaDiversity.xlsx")))

# ----------------------------------------
# Top contributing species based on absolute difference
# ----------------------------------------
species_diff <- abs(species_matrix["PatientReads", ] - species_matrix["ControlReads", ])
top_beta <- sort(species_diff, decreasing = TRUE)[1:10]
top_beta_df <- data.frame(Species = names(top_beta), AbsDiff = top_beta)

# Plot top contributing species
beta_plot <- ggplot(top_beta_df, aes(x = AbsDiff, y = reorder(Species, AbsDiff))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Top Contributing Species to Bray-Curtis Distance",
       x = "Absolute Difference",
       y = "Species") +
  theme_minimal()

# Save plot
beta_file <- file.path(outdir, paste0(patient_name, "_BetaDiversity.png"))
ggsave(beta_file, beta_plot, width = 8, height = 6, dpi = 300)

# =========================================================
# 10. Alpha Diversity
# =========================================================
alpha_div <- data.frame(
  Sample=rownames(species_matrix),
  Observed=apply(species_matrix,1,function(x) sum(x>0)),
  Shannon=vegan::diversity(species_matrix, index="shannon"),
  Simpson=vegan::diversity(species_matrix, index="simpson")
)
write_xlsx(alpha_div, file.path(outdir,paste0(patient_name,"_AlphaDiversity.xlsx")))

# Visualize alpha diversity
alpha_long <- alpha_div %>% pivot_longer(cols=-Sample, names_to="Metric", values_to="Value")
alpha_plot <- ggplot(alpha_long, aes(x=Sample, y=Value, fill=Sample)) +
  geom_bar(stat="identity", position="dodge") + facet_wrap(~Metric, scales="free_y") +
  labs(title="Alpha Diversity Metrics", x="Sample", y="Value") + theme_minimal() +
  theme(legend.position="none")
ggsave(file.path(outdir, paste0(patient_name,"_AlphaDiversity.png")), alpha_plot, width=10, height=6, dpi=300)
