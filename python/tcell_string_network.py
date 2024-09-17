
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load the T-cell interaction network file
tcell_interactions_path = './../data/Tcell_interactions.tsv'
tcell_interactions = pd.read_csv(tcell_interactions_path, sep="\t")

# Create a directed graph using the interaction data
G = nx.DiGraph()

# Add edges with combined score as edge weight
for index, row in tcell_interactions.iterrows():
    G.add_edge(row['#node1'], row['node2'], weight=row['combined_score'])

# Visualize the graph
plt.figure(figsize=(12, 8))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color="lightblue", edge_color="gray", node_size=500, font_size=10)
plt.title("T-cell STRING Network")
plt.show()
