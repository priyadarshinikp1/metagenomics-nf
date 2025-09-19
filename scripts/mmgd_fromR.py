#!/usr/bin/env python3
# py mmgd_fromR.py p975_Microbiome_Filtered_RelativeAbundance.xlsx "Gut Microbe-Microbial metabolite.csv" "Microbial metabolite-Host Gene.csv" "Gut Microbe-Host Gene.csv"

import pandas as pd
import numpy as np
import sys
import openpyxl
import os

def load_table(file):
    """Try reading CSV or Excel automatically with correct encoding"""
    try:
        return pd.read_csv(file, encoding="utf-8")
    except UnicodeDecodeError:
        try:
            return pd.read_csv(file, encoding="latin1")
        except Exception:
            return pd.read_excel(file)

if len(sys.argv) != 5:
    print(f"Usage: {sys.argv[0]} <R_diff_abundance.xlsx> "
          "<gutmgene_microbe_metabolite.csv> <gutmgene_metabolite_gene.csv> "
          "<gutmgene_microbe_gene.csv>")
    sys.exit(1)

r_diff_file, mm_file, mg_file, mhg_file = sys.argv[1:]

# ---------------- Output prefix ----------------
full_prefix = os.path.splitext(os.path.basename(r_diff_file))[0]
prefix = full_prefix.split("_")[0]  # keep only first token, e.g. "p975"


# ---------------- Differential Abundance from R ----------------
diff = pd.read_excel(r_diff_file)

# Standardize column names coming from your R export
df = diff.rename(columns={
    "Name": "Matched_Microbe",
    "Reads_Clade": "patient_RA",
    "log2FoldChange": "Log2FC"
})

# Add helper fields
df["Rank"] = "S"
df["clean_taxon"] = df["Matched_Microbe"].str.lower()

# ---------------- GutMGene ----------------
mm = load_table(mm_file)   # microbe-metabolite
mg = load_table(mg_file)   # metabolite-gene
mhg = load_table(mhg_file) # microbe-gene

mm["clean_microbe"] = mm["Gut Microbiota"].str.lower()
mhg["clean_microbe"] = mhg["Gut Microbiota"].str.lower()

# join microbe-metabolite
merged = pd.merge(df, mm, left_on="clean_taxon", right_on="clean_microbe", how="left")

# join metabolite-gene
merged = pd.merge(merged, mg, on="Metabolite", how="left", suffixes=("", "_mg"))

# join microbe-host gene
merged = pd.merge(merged, mhg, left_on="clean_taxon", right_on="clean_microbe", how="left", suffixes=("", "_mhg"))

# ---------------- Final table ----------------
final = merged[[
    "Matched_Microbe", "Rank", "Metabolite", "Gene", "Alteration",
    "patient_RA", "Log2FC", "Condition"
]].drop_duplicates()

final = final.rename(columns={"Condition": "Disease"})

# keep only matched rows (drop empty mappings)
final = final[(final["Metabolite"].notna()) &
              (final["Gene"].notna()) &
              (final["Disease"].notna())]

# filter out healthy conditions
final = final[final["Disease"].str.lower() != "healthy"]

final["Match_Status"] = "matched"

# ---------------- Save CSV ----------------
csv_output = f"{prefix}_microbe_metabolite_gene_disease_network.csv"
final.to_csv(csv_output, index=False)
print(f"[INFO] Final matched results written to {csv_output}")

# ---------------- Network visualization ----------------
import networkx as nx
from pyvis.network import Network

top_interactions = final.copy()
top_interactions['abs_Log2FC'] = top_interactions['Log2FC'].abs()
top_interactions = top_interactions.sort_values(by='abs_Log2FC', ascending=False).head(100)

net = Network(height='750px', width='100%', notebook=False, directed=False)
net.set_options("""
var options = {
  "physics": {
    "barnesHut": {
      "gravitationalConstant": -3000,
      "springLength": 150,
      "springConstant": 0.05,
      "damping": 0.75
    },
    "minVelocity": 0.1,
    "stabilization": {
      "enabled": true,
      "iterations": 1000,
      "updateInterval": 50,
      "onlyDynamicEdges": false,
      "fit": true
    }
  }
}
""")

# Add nodes
microbes = top_interactions['Matched_Microbe'].unique()
metabolites = top_interactions['Metabolite'].unique()
genes = top_interactions['Gene'].unique()
diseases = top_interactions['Disease'].unique()

for m in microbes:
    net.add_node(m, label=m, color='lightblue', title='Microbe')
for met in metabolites:
    net.add_node(met, label=met, color='orange', title='Metabolite')
for g in genes:
    net.add_node(g, label=g, color='green', title='Gene')
for d in diseases:
    net.add_node(d, label=d, color='red', title='Disease')

# Add edges
for _, row in top_interactions.iterrows():
    net.add_edge(row['Matched_Microbe'], row['Metabolite'])
    net.add_edge(row['Metabolite'], row['Gene'])
    net.add_edge(row['Matched_Microbe'], row['Gene'])
    net.add_edge(row['Gene'], row['Disease'])

# ---------------- Save interactive HTML ----------------
html_output = f"{prefix}_top_microbe_metabolite_gene_disease_network.html"
net.write_html(html_output)
print(f"[INFO] Network visualization saved as {html_output}")

# ---------------- Save static PNG using matplotlib ----------------
import matplotlib.pyplot as plt

G = nx.Graph()
for _, row in top_interactions.iterrows():
    G.add_edge(row['Matched_Microbe'], row['Metabolite'])
    G.add_edge(row['Metabolite'], row['Gene'])
    G.add_edge(row['Matched_Microbe'], row['Gene'])
    G.add_edge(row['Gene'], row['Disease'])

plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, k=0.3, seed=42)  # seed for reproducibility

# Draw nodes by category
nx.draw_networkx_nodes(G, pos, nodelist=microbes, node_color="lightblue", node_size=500, label="Microbes")
nx.draw_networkx_nodes(G, pos, nodelist=metabolites, node_color="orange", node_size=500, label="Metabolites")
nx.draw_networkx_nodes(G, pos, nodelist=genes, node_color="green", node_size=500, label="Genes")
nx.draw_networkx_nodes(G, pos, nodelist=diseases, node_color="red", node_size=500, label="Diseases")

# Draw edges + labels
nx.draw_networkx_edges(G, pos, alpha=0.3)
nx.draw_networkx_labels(G, pos, font_size=8)

plt.legend(scatterpoints=1)
plt.axis("off")

png_output = f"{prefix}_top_microbe_metabolite_gene_disease_network.png"
plt.savefig(png_output, dpi=300, bbox_inches="tight")
plt.close()
print(f"[INFO] Static PNG network saved as {png_output}")
