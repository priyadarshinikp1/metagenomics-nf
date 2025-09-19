# ==========================================
# Complete Microbe-Metabolite-Gene-Disease Analysis
# With Automatic CSV, Plot, Interactive HTML Network, Sankey, and Report
# ==========================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from gprofiler import GProfiler
import plotly.graph_objects as go
import os
import plotly.io as pio
from pyvis.network import Network
from jinja2 import Template

sns.set(style="whitegrid")

# =========================
# 1. Load Data
# =========================
import glob, sys

csv_files = glob.glob("*_microbe_metabolite_gene_disease_network.csv")

if not csv_files:
    print("[ERROR] No matching CSV file found (expected '*_microbe_metabolite_gene_disease_network.csv').")
    sys.exit(1)

# If multiple exist, just take the first
input_file = csv_files[0]
print(f"[INFO] Using input file: {input_file}")

df = pd.read_csv(input_file)

# derive prefix from file name (before "_microbe_metabolite_gene_disease_network")
prefix = os.path.splitext(os.path.basename(input_file))[0].split("_")[0]


# =========================
# 2. Data Cleaning
# =========================
df['Gene'] = df['Gene'].str.upper()
df['Metabolite'].fillna('Unknown', inplace=True)
df['Matched_Microbe'].fillna('Unknown', inplace=True)
df['Log2FC'] = pd.to_numeric(df['Log2FC'], errors='coerce')
df.drop_duplicates(inplace=True)

# =========================
# 3. Create Output Folder
# =========================
output_folder = "analysis_outputs"
os.makedirs(output_folder, exist_ok=True)

# =========================
# 4. Summary Tables
# =========================
top_microbes = df['Matched_Microbe'].value_counts().reset_index()
top_microbes.columns = ['Microbe','Num_Interactions']
top_microbes.to_csv(f"{output_folder}/{prefix}_top_microbes.csv", index=False)

top_metabolites = df['Metabolite'].value_counts().reset_index()
top_metabolites.columns = ['Metabolite','Num_Interactions']
top_metabolites.to_csv(f"{output_folder}/{prefix}_top_metabolites.csv", index=False)

alteration_counts = df.groupby(['Metabolite','Alteration']).size().unstack(fill_value=0)
alteration_counts.to_csv(f"{output_folder}/{prefix}_metabolite_alteration_counts.csv")

disease_microbe = df.groupby(['Disease','Matched_Microbe'])['Gene'].nunique().reset_index()
disease_microbe.rename(columns={'Gene':'Num_Genes_Affected'}, inplace=True)
top_microbes_per_disease = disease_microbe.sort_values(['Disease','Num_Genes_Affected'], ascending=[True, False])
top_microbes_per_disease.to_csv(f"{output_folder}/{prefix}_top_microbes_per_disease.csv", index=False)

disease_metabolite = df.groupby(['Disease','Metabolite'])['Gene'].nunique().reset_index()
disease_metabolite.rename(columns={'Gene':'Num_Genes_Affected'}, inplace=True)
top_metabolites_per_disease = disease_metabolite.sort_values(['Disease','Num_Genes_Affected'], ascending=[True, False])
top_metabolites_per_disease.to_csv(f"{output_folder}/{prefix}_top_metabolites_per_disease.csv", index=False)

disease_metabolite_fc = df.groupby(['Disease','Metabolite'])['Log2FC'].mean().reset_index()
disease_metabolite_fc.to_csv(f"{output_folder}/{prefix}_top_metabolites_per_disease_log2fc.csv", index=False)

# =========================
# 5. Visualizations
# =========================
plt.figure(figsize=(10,6))
sns.barplot(x=top_microbes['Num_Interactions'], y=top_microbes['Microbe'], palette="viridis")
plt.title("Top 10 Microbes by Number of Interactions")
plt.xlabel("Number of Interactions")
plt.ylabel("Microbe")
plt.tight_layout()
plt.savefig(f"{output_folder}/{prefix}_top_microbes_barplot.png")
plt.close()

plt.figure(figsize=(10,6))
sns.barplot(x=top_metabolites['Num_Interactions'], y=top_metabolites['Metabolite'], palette="magma")
plt.title("Top 10 Metabolites by Number of Gene Interactions")
plt.xlabel("Number of Genes Affected")
plt.ylabel("Metabolite")
plt.tight_layout()
plt.savefig(f"{output_folder}/{prefix}_top_metabolites_barplot.png")
plt.close()

plt.figure(figsize=(12,6))
sns.heatmap(alteration_counts, annot=True, fmt="d", cmap="coolwarm")
plt.title("Metabolite Alteration Counts (Activation vs Inhibition)")
plt.tight_layout()
plt.savefig(f"{output_folder}/{prefix}_metabolite_alteration_heatmap.png")
plt.close()

disease_gene_matrix = df.pivot_table(index='Gene', columns='Disease', values='Log2FC', aggfunc='mean', fill_value=0)
plt.figure(figsize=(12,10))
sns.heatmap(disease_gene_matrix, cmap="RdBu_r", center=0)
plt.title("Gene Log2FC Across Diseases")
plt.tight_layout()
plt.savefig(f"{output_folder}/{prefix}_disease_gene_log2fc_heatmap.png")
plt.close()

# =========================
# 6. Pathway Enrichment Analysis
# =========================
gp = GProfiler(return_dataframe=True)
disease_list = df['Disease'].unique()
pathway_results_all = pd.DataFrame()

for disease in disease_list:
    genes = df[df['Disease']==disease]['Gene'].unique()
    if len(genes) < 3:
        continue
    res = gp.profile(organism='mmusculus', query=list(genes), sources=['GO:BP','KEGG','REAC'])
    res['Disease'] = disease
    pathway_results_all = pd.concat([pathway_results_all, res], ignore_index=True)

pathway_results_all.to_csv(f"{output_folder}/{prefix}_pathway_enrichment_all_diseases.csv", index=False)

for disease in disease_list:
    res = pathway_results_all[pathway_results_all['Disease']==disease].head(10)
    if res.empty:
        continue
    plt.figure(figsize=(10,6))
    sns.barplot(x=-np.log10(res['p_value']), y=res['name'], palette="viridis")
    plt.xlabel("-log10(p-value)")
    plt.ylabel("Pathway")
    plt.title(f"Top Pathways for {disease}")
    plt.tight_layout()
    plt.savefig(f"{output_folder}/{prefix}_top_pathways_{disease.replace(' ','_')}.png")
    plt.close()

# =========================
# 7. Network Visualization (PNG + HTML)
# =========================
G = nx.DiGraph()
for _, row in df.iterrows():
    G.add_edge(row['Matched_Microbe'], row['Metabolite'], type='microbe_metabolite')
    G.add_edge(row['Metabolite'], row['Gene'], type=row['Alteration'])
    G.add_edge(row['Gene'], row['Disease'], type='gene_disease')

plt.figure(figsize=(15,15))
pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
edge_colors = [
    'green' if G[u][v]['type'].lower() == 'activation' 
    else 'red' if G[u][v]['type'].lower() == 'inhibition' 
    else 'gray' 
    for u, v in G.edges()
]
nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=500, edge_color=edge_colors, font_size=8)
plt.title("Microbe-Metabolite-Gene-Disease Network")
plt.tight_layout()
png_network = f"{output_folder}/{prefix}_network_graph.png"
plt.savefig(png_network, dpi=300, bbox_inches="tight")
plt.close()
print(f"[INFO] Static PNG network saved as {png_network}")

net = Network(height="750px", width="100%", directed=True, notebook=False)
for node in G.nodes():
    if node in df['Matched_Microbe'].values:
        net.add_node(node, label=node, color="skyblue", title=f"Microbe: {node}")
    elif node in df['Metabolite'].values:
        net.add_node(node, label=node, color="orange", title=f"Metabolite: {node}")
    elif node in df['Gene'].values:
        net.add_node(node, label=node, color="green", title=f"Gene: {node}")
    elif node in df['Disease'].values:
        net.add_node(node, label=node, color="red", title=f"Disease: {node}")
    else:
        net.add_node(node, label=node, color="gray")
for u, v, d in G.edges(data=True):
    edge_type = d["type"].lower()
    color = "green" if edge_type == "activation" else "red" if edge_type == "inhibition" else "gray"
    net.add_edge(u, v, color=color, title=edge_type)
html_network = f"{output_folder}/{prefix}_network_graph.html"
net.write_html(html_network)
print(f"[INFO] Interactive HTML network saved as {html_network}")

# =========================
# 8. Sankey Diagram
# =========================
microbes = list(df['Matched_Microbe'].unique())
metabolites = list(df['Metabolite'].unique())
genes = list(df['Gene'].unique())
diseases = list(df['Disease'].unique())
nodes = microbes + metabolites + genes + diseases
node_indices = {node: i for i, node in enumerate(nodes)}

sources, targets, values, colors = [], [], [], []
for _, row in df.iterrows():
    sources.append(node_indices[row['Matched_Microbe']]); targets.append(node_indices[row['Metabolite']]); values.append(1); colors.append('#add8e6')
    sources.append(node_indices[row['Metabolite']]); targets.append(node_indices[row['Gene']]); values.append(1); colors.append('#90ee90' if row['Alteration'].lower()=='activation' else '#f08080')
    sources.append(node_indices[row['Gene']]); targets.append(node_indices[row['Disease']]); values.append(1); colors.append('#90ee90' if row['Alteration'].lower()=='activation' else '#f08080')

sankey_fig = go.Figure(go.Sankey(
    node=dict(pad=20, thickness=25, label=nodes, color="skyblue", line=dict(color="black", width=0.5), hovertemplate='%{label}<extra></extra>'),
    link=dict(source=sources, target=targets, value=values, color=colors)
))
sankey_fig.update_layout(title_text="Microbe-Metabolite-Gene-Disease Sankey Diagram", font_size=12)
sankey_div = pio.to_html(sankey_fig, full_html=False)
with open(f"{output_folder}/{prefix}_sankey_diagram_colored.html", "w", encoding="utf-8") as f:
    f.write(sankey_div)
print("Interactive colored Sankey diagram saved!")

# =========================
# 9. HTML Report
# =========================
top_microbes = pd.read_csv(f"{output_folder}/{prefix}_top_microbes.csv")
top_metabolites = pd.read_csv(f"{output_folder}/{prefix}_top_metabolites.csv")
top_microbes_disease = pd.read_csv(f"{output_folder}/{prefix}_top_microbes_per_disease.csv")
top_metabolites_disease = pd.read_csv(f"{output_folder}/{prefix}_top_metabolites_per_disease.csv")
pathways = pd.read_csv(f"{output_folder}/{prefix}_pathway_enrichment_all_diseases.csv", encoding="utf-8")
pathway_imgs = [f for f in os.listdir(output_folder) if f.startswith(f"{prefix}_top_pathways_") and f.endswith(".png")]

plots = {
    "Top Microbes": f"{prefix}_top_microbes_barplot.png",
    "Top Metabolites": f"{prefix}_top_metabolites_barplot.png",
    "Metabolite Alteration Heatmap": f"{prefix}_metabolite_alteration_heatmap.png",
    "Gene-Disease Log2FC Heatmap": f"{prefix}_disease_gene_log2fc_heatmap.png",
}

html_template = """
<html>
<head>
    <title>Microbe-Metabolite-Gene-Disease Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 0; padding: 0; background: #f9f9f9; }
        .container { max-width: 1000px; margin: auto; padding: 20px; background: white; }
        h1,h2,h3 { color: #2C3E50; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        img { max-width: 100%; height: auto; display: block; margin: 10px auto; }
        details { margin-bottom: 30px; }
        summary { font-size: 1.2em; font-weight: bold; cursor: pointer; }
    </style>
</head>
<body>
<div class="container">
    <h1>Microbe-Metabolite-Gene-Disease Analysis Report</h1>

    <details open><summary>Top Microbes</summary>{{ top_microbes_table }}<img src="{{ plots['Top Microbes'] }}"></details>
    <details open><summary>Top Metabolites</summary>{{ top_metabolites_table }}<img src="{{ plots['Top Metabolites'] }}"></details>
    <details><summary>Metabolite Alteration Counts</summary><img src="{{ plots['Metabolite Alteration Heatmap'] }}"></details>
    <details><summary>Gene-Disease Log2FC Heatmap</summary><img src="{{ plots['Gene-Disease Log2FC Heatmap'] }}"></details>
    <details><summary>Top Microbes per Disease</summary>{{ top_microbes_disease_table }}</details>
    <details><summary>Top Metabolites per Disease</summary>{{ top_metabolites_disease_table }}</details>

    <details open><summary>Top Pathway Enrichments (Table)</summary>
        {{ pathways_table }}
    </details>

    <details open><summary>Top Pathway Enrichments (Plots)</summary>
        {% for img in pathway_imgs %}
          <div class="pathway-section">
            <img src="{{ img }}">
          </div>
        {% endfor %}
    </details>

    <details open><summary>Interactive Colored Sankey Diagram</summary>{{ sankey_div | safe }}</details>
</div>
</body>
</html>
"""

template = Template(html_template)
html_content = template.render(
    top_microbes_table=top_microbes.head(10).to_html(index=False),
    top_metabolites_table=top_metabolites.head(10).to_html(index=False),
    top_microbes_disease_table=top_microbes_disease.head(10).to_html(index=False),
    top_metabolites_disease_table=top_metabolites_disease.head(10).to_html(index=False),
    pathways_table=pathways.head(5).to_html(index=False),   # <-- added table
    plots=plots,
    pathway_imgs=pathway_imgs,
    sankey_div=sankey_div
)

report_file = f"{output_folder}/{prefix}_microbe_metabolite_gene_disease_report.html"
with open(report_file, "w", encoding="utf-8") as f:
    f.write(html_content)
print(f"Final report with Sankey and pathway tables saved to: {report_file}")
