"""Pathway enrichment analysis on the real SARS-CoV-2 host-interactome (Gordon et al. 2020),
using Enrichr (GO Biological Process, KEGG, Reactome) via gseapy.
"""
import pandas as pd
import gseapy as gp

LIBRARIES = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
MIN_GENES = 3


def enrich(genes, label):
    try:
        res = gp.enrichr(gene_list=genes, gene_sets=LIBRARIES, organism="human", outdir=None)
        df = res.results
        df = df[df["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value")
        df.insert(0, "query", label)
        return df
    except Exception as ex:
        print(f"  [skip] {label}: {ex}")
        return pd.DataFrame()


def main():
    df = pd.read_csv("data/interactome.csv")

    # 1. Overall pathway landscape across the whole interactome
    all_genes = sorted(df["prey_gene"].unique().tolist())
    print(f"Overall enrichment across {len(all_genes)} host proteins...")
    overall = enrich(all_genes, "ALL_BAITS")
    overall.to_csv("data/pathways_overall.csv", index=False)
    print(f"  -> {len(overall)} significant terms (FDR<0.05)")

    # 2. Per-bait enrichment (which pathways each viral protein's targets fall into)
    per_bait_frames = []
    for bait, sub in df.groupby("bait"):
        genes = sub["prey_gene"].dropna().unique().tolist()
        if len(genes) < MIN_GENES:
            continue
        print(f"Enrichment for {bait} ({len(genes)} preys)...")
        res = enrich(genes, bait)
        if not res.empty:
            per_bait_frames.append(res)

    per_bait = pd.concat(per_bait_frames, ignore_index=True) if per_bait_frames else pd.DataFrame()
    per_bait.to_csv("data/pathways_per_bait.csv", index=False)
    print(f"Done. {len(per_bait)} significant bait-pathway hits saved to data/pathways_per_bait.csv")


if __name__ == "__main__":
    main()
