"""Visualize two predicted-missing-component examples as small interaction graphs:
  1. Nsp13 -> PLK1 (strongest single guilt-by-association call, centrosome pathway)
  2. Nsp4 -> TOMM70, cross-validated because TOMM70 is a REAL confirmed Orf9b interactor
     in the paper's own data -- independent convergence on the same mitochondrial receptor.
"""
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("data/interactome.csv")
pred = pd.read_csv("data/predicted_missing_components.csv")


def plot_case(bait, candidate, extra_bait=None, extra_known_edge=None, fname="pred.png", title=""):
    row = pred[(pred["bait"] == bait) & (pred["candidate_gene"] == candidate)].iloc[0]
    connected = row["connected_to"].split(";")

    G = nx.Graph()
    G.add_node(bait, kind="bait")
    for g in connected:
        G.add_node(g, kind="known_prey")
        G.add_edge(bait, g, kind="known")
    G.add_node(candidate, kind="predicted")
    for g in connected:
        G.add_edge(candidate, g, kind="predicted")

    if extra_bait and extra_known_edge:
        G.add_node(extra_bait, kind="bait")
        G.add_edge(extra_bait, candidate, kind="known")

    pos = nx.spring_layout(G, seed=7, k=0.9)
    colors, edge_colors, styles = [], [], []
    for n, d in G.nodes(data=True):
        colors.append({"bait": "#c0392b", "known_prey": "#3b6fa0", "predicted": "#e08e2b"}[d["kind"]])
    for u, v, d in G.edges(data=True):
        if d["kind"] == "known":
            edge_colors.append("#666666"); styles.append("solid")
        else:
            edge_colors.append("#e08e2b"); styles.append("dashed")

    fig, ax = plt.subplots(figsize=(7, 6))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors, style=styles, width=1.8)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=colors, node_size=900)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)
    ax.set_title(title, fontsize=10)
    ax.axis("off")
    from matplotlib.lines import Line2D
    legend = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#c0392b', markersize=10, label='Viral bait'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#3b6fa0', markersize=10, label='Known host interactor'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#e08e2b', markersize=10, label='Predicted missing component'),
    ]
    ax.legend(handles=legend, loc="upper right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print("saved", fname)


plot_case("Nsp13", "PLK1", fname="output/predicted_nsp13_plk1.png",
          title="Nsp13 centrosome network: PLK1 predicted as missing component\n"
                "(guilt-by-association via STRING, 15 edges to known Nsp13 preys)")

plot_case("Nsp4", "TOMM70", extra_bait="Orf9b", extra_known_edge=True,
          fname="output/predicted_nsp4_tomm70.png",
          title="TOMM70: predicted new Nsp4 interactor, cross-validated\n"
                "TOMM70 is Orf9b's REAL confirmed interactor in the paper's own data")
