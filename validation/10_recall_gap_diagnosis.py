"""Next-step 3: diagnose the bait-dependent recall gap.

Test 3's internal holdout recovered held-out interactions for some baits (N: 46.7%,
Nsp13: 27.5%) and none at all for others (M, Nsp7, Nsp12, Orf9b: 0%). Working
hypothesis: the STRING guilt-by-association method only works when a bait's known
interactors already form a tight, densely-connected cluster in STRING's evidence graph
-- so add_nodes has a coherent neighborhood to extend -- and that this pre-existing
density is systematically thinner for membrane/trafficking baits than for well-studied
complexes.

This script quantifies that directly. For each bait tested in the holdout, it pulls the
STRING sub-network AMONG that bait's known preys only (no add_nodes), measures how dense
that intra-set connectivity is, and -- crucially -- decomposes the density by evidence
channel, so we can see whether the baits that "work" are held together by hard
experimental evidence (escore) or by literature-biased channels (tscore text-mining,
dscore curated database) that inflate well-studied biology.

Then it correlates each density metric against the bait's mean holdout recall (test 3).
If recall tracks intra-set density, that explains the gap; the channel breakdown says
whether to trust it.

STRING channels: nscore (neighborhood), fscore (fusion), pscore (phylo co-occurrence),
ascore (coexpression), escore (experiments), dscore (curated database), tscore (text-mining).
"""
import time
import math
import urllib.request
import urllib.parse
from io import StringIO
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STRING_URL = "https://string-db.org/api/tsv/network"
CHANNELS = ["escore", "dscore", "tscore", "ascore", "nscore", "fscore", "pscore"]
# escore = hard experimental evidence; tscore/dscore = literature/curation (can be inflated
# for heavily-studied complexes); ascore = coexpression.


def string_intra_network(genes):
    """STRING edges AMONG the given genes only (no expansion)."""
    ids = "%0d".join(genes)
    url = STRING_URL + "?" + urllib.parse.urlencode(
        {"identifiers": ids, "species": "9606"}, safe="%0d")
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()
    if not text.strip():
        return pd.DataFrame()
    return pd.read_csv(StringIO(text), sep="\t")


def main():
    gordon = pd.read_csv("data/interactome.csv")

    # per-bait mean real recall from the internal holdout (test 3)
    holdout = pd.read_csv("data/validation/holdout_precision_recall.csv")
    recall = holdout[holdout["condition"] == "real"].groupby("bait")["recall"].mean()

    rows = []
    for bait in recall.index:
        genes = sorted(gordon.loc[gordon["bait"] == bait, "prey_gene"].tolist())
        n = len(genes)
        n_pairs = math.comb(n, 2)
        print(f"{bait}: {n} preys, {n_pairs} possible intra-set pairs -> querying STRING...")
        edges = string_intra_network(genes)
        n_edges = len(edges)

        # density-weighted channel means: sum(channel score over present edges) / all possible
        # pairs. Captures BOTH how many preys are linked and how strong each evidence type is.
        rec = {
            "bait": bait,
            "n_preys": n,
            "n_intra_edges": n_edges,
            "edge_density": n_edges / n_pairs if n_pairs else float("nan"),
            "mean_combined_over_pairs": (edges["score"].sum() / n_pairs) if (n_pairs and n_edges) else 0.0,
            "holdout_recall": recall[bait],
        }
        for ch in CHANNELS:
            rec[f"{ch}_density"] = (edges[ch].sum() / n_pairs) if (n_pairs and n_edges) else 0.0
        rows.append(rec)
        time.sleep(1)

    df = pd.DataFrame(rows).sort_values("holdout_recall", ascending=False)
    df.to_csv("data/validation/recall_gap_diagnosis.csv", index=False)

    print("\n=== Per-bait intra-set STRING density vs holdout recall ===")
    show = ["bait", "n_preys", "edge_density", "escore_density", "tscore_density",
            "dscore_density", "ascore_density", "holdout_recall"]
    print(df[show].round(3).to_string(index=False))

    # correlations (Spearman -- monotonic, robust to the small n and skew)
    print("\n=== Spearman correlation of each density metric with holdout recall ===")
    metrics = ["edge_density", "mean_combined_over_pairs"] + [f"{ch}_density" for ch in CHANNELS]
    corrs = {m: df[m].corr(df["holdout_recall"], method="spearman") for m in metrics}
    for m, c in sorted(corrs.items(), key=lambda x: -(x[1] if not np.isnan(x[1]) else -9)):
        print(f"  {m:28s} rho = {c:+.3f}")

    # figure: intra-set edge density vs recall, and the experiments-vs-textmining split
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4))

    ax1.scatter(df["edge_density"], df["holdout_recall"], s=90, color="#3b6fa0", zorder=3)
    for _, r in df.iterrows():
        ax1.annotate(r["bait"], (r["edge_density"], r["holdout_recall"]),
                     xytext=(5, 4), textcoords="offset points", fontsize=9)
    ax1.set_xlabel("Intra-prey-set STRING edge density\n(fraction of possible prey-prey pairs connected)")
    ax1.set_ylabel("Mean holdout recall (test 3)")
    ax1.set_title("Does the method work only where the bait's\nbiology is already a tight STRING cluster?")
    ax1.grid(True, alpha=0.25)

    # experiments-channel vs text-mining-channel density, point color = recall
    sc = ax2.scatter(df["escore_density"], df["tscore_density"], c=df["holdout_recall"],
                     cmap="viridis", s=110, zorder=3, edgecolor="white", linewidth=0.5)
    for _, r in df.iterrows():
        ax2.annotate(r["bait"], (r["escore_density"], r["tscore_density"]),
                     xytext=(5, 4), textcoords="offset points", fontsize=9)
    lim = max(df["escore_density"].max(), df["tscore_density"].max()) * 1.1
    ax2.plot([0, lim], [0, lim], ls="--", color="#999", lw=1)
    ax2.set_xlabel("Experiments-channel density (escore)\nhard experimental evidence")
    ax2.set_ylabel("Text-mining-channel density (tscore)\nliterature-biased")
    ax2.set_title("Is the connectivity driven by hard evidence\nor by literature co-mention?")
    ax2.grid(True, alpha=0.25)
    cbar = fig.colorbar(sc, ax=ax2)
    cbar.set_label("holdout recall")

    fig.tight_layout()
    fig.savefig("output/validation/recall_gap_diagnosis.png", dpi=150)
    print("\nsaved output/validation/recall_gap_diagnosis.png")
    print("saved data/validation/recall_gap_diagnosis.csv")


if __name__ == "__main__":
    main()
