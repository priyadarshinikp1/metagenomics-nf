import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import os
import glob
import sys

# ---------------- Load input automatically ----------------
csv_files = glob.glob("*_Microbe_Gene_RankSimilarity.csv")

if not csv_files:
    print("[ERROR] No matching CSV file found (expected '*_Microbe_Gene_RankSimilarity.csv').")
    sys.exit(1)

# If multiple exist, just take the first
input_file = csv_files[0]
print(f"[INFO] Using input file: {input_file}")

df = pd.read_csv(input_file)

# Prefix for outputs (before "_Microbe_Gene_RankSimilarity")
prefix = os.path.splitext(os.path.basename(input_file))[0]
prefix = prefix.replace("-", "_").split("_")[0]


# ---------------- Filter data ----------------
# Top 20 microbes by MicrobeRank
top_microbes = (
    df.sort_values('MicrobeRank')
      .drop_duplicates('Microbe_clean')
      .head(20)['Microbe_clean']
      .tolist()
)

# Top 100 strongest edges among those microbes
filtered_df = (
    df[df['Microbe_clean'].isin(top_microbes)]
      .sort_values('RankSimilarity', ascending=False)
      .head(100)
)

# ---------------- PyVis interactive network ----------------
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
microbes = filtered_df['Microbe_clean'].unique()
genes = filtered_df['GeneID'].unique()

for m in microbes:
    microbe_rank = df.loc[df['Microbe_clean'] == m, 'MicrobeRank'].values[0]
    net.add_node(m, label=m, color='lightblue',
                 title=f"Microbe: {m}\nRank: {microbe_rank}", size=25)

for g in genes:
    net.add_node(g, label=g, color='green',
                 title=f"Gene: {g}", size=15)

# Add edges
for _, row in filtered_df.iterrows():
    net.add_edge(row['Microbe_clean'], row['GeneID'],
                 value=row['RankSimilarity'] * 10, color='blue')

# ---------------- Save interactive HTML ----------------
html_output = f"{prefix}_Microbe_Gene_rank_network.html"
net.write_html(html_output)
print(f"[INFO] Interactive network saved as {html_output}")

# ---------------- Save static PNG (matplotlib) ----------------
G = nx.Graph()
for _, row in filtered_df.iterrows():
    G.add_edge(row['Microbe_clean'], row['GeneID'], weight=row['RankSimilarity'])

plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, k=0.3, seed=42)  # reproducible layout

# Draw nodes
nx.draw_networkx_nodes(G, pos, nodelist=microbes,
                       node_color="lightblue", node_size=600, label="Microbes")
nx.draw_networkx_nodes(G, pos, nodelist=genes,
                       node_color="green", node_size=400, label="Genes")

# Draw edges weighted by similarity
edges = nx.draw_networkx_edges(
    G, pos,
    width=[d['weight'] * 5 for (_, _, d) in G.edges(data=True)],
    alpha=0.4
)

# Labels
nx.draw_networkx_labels(G, pos, font_size=8)

plt.legend(scatterpoints=1)
plt.axis("off")

png_output = f"{prefix}_Microbe_Gene_rank_network.png"
plt.savefig(png_output, dpi=300, bbox_inches="tight")
plt.close()
print(f"[INFO] Static PNG network saved as {png_output}")
