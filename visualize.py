"""Recreate a Fig.2b-style bait x pathway enrichment heatmap, plus an overall top-pathways bar chart,
from real enrichment results computed on the Gordon et al. 2020 SARS-CoV-2 interactome."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BAIT_ORDER = ["Nsp1", "Nsp2", "Nsp4", "Nsp5", "Nsp5 C145A", "Nsp6", "Nsp7", "Nsp8", "Nsp9",
              "Nsp10", "Nsp11", "Nsp12", "Nsp13", "Nsp14", "Nsp15", "Spike", "E", "M", "N",
              "Orf3a", "Orf3b", "Orf6", "Orf7a", "Orf8", "Orf9b", "Orf9c", "Orf10"]


def plot_overall_bar():
    df = pd.read_csv("data/pathways_overall.csv").head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 6))
    vals = -np.log10(df["Adjusted P-value"])
    ax.barh(df["Term"].str.slice(0, 55), vals, color="#3b6fa0")
    ax.set_xlabel("-log10(adjusted P-value)")
    ax.set_title("Top enriched pathways across all 332 SARS-CoV-2 host-interacting proteins\n"
                  "(real data, Gordon et al. 2020 interactome + Enrichr)")
    fig.tight_layout()
    fig.savefig("output/overall_top_pathways.png", dpi=150)
    print("saved output/overall_top_pathways.png")


def plot_bait_pathway_heatmap(top_n_terms=25):
    df = pd.read_csv("data/pathways_per_bait.csv")
    df["neglogq"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))

    # pick the globally most significant terms to keep the heatmap readable
    best_per_term = df.groupby("Term")["neglogq"].max().sort_values(ascending=False)
    top_terms = best_per_term.head(top_n_terms).index.tolist()

    baits_present = [b for b in BAIT_ORDER if b in df["query"].unique()]
    mat = pd.DataFrame(0.0, index=top_terms, columns=baits_present)
    for _, row in df[df["Term"].isin(top_terms)].iterrows():
        mat.loc[row["Term"], row["query"]] = max(mat.loc[row["Term"], row["query"]], row["neglogq"])

    fig, ax = plt.subplots(figsize=(13, 9))
    im = ax.imshow(mat.values, cmap="Blues", aspect="auto", vmin=0, vmax=mat.values.max())
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=90)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels([t[:60] for t in mat.index], fontsize=8)
    cbar = fig.colorbar(im, ax=ax, shrink=0.6)
    cbar.set_label("-log10(FDR)")
    ax.set_title("Pathway enrichment per SARS-CoV-2 viral protein (real data)\n"
                  "reproducing the Fig. 2b logic of Gordon et al. 2020")
    fig.tight_layout()
    fig.savefig("output/bait_pathway_heatmap.png", dpi=150)
    print("saved output/bait_pathway_heatmap.png")


if __name__ == "__main__":
    plot_overall_bar()
    plot_bait_pathway_heatmap()
