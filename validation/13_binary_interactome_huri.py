"""Next-step (reviewer request): test Nsp13->PLK1 against a NON-curated binary interactome.

Tests 8, 9, and 12 all leaned on literature-curated data (STRING, IntAct, published AP-MS).
Even the positive cross-method result (test 12, IntAct) carries the caveat that IntAct is
literature-curated, so heavily-studied centrosome biology is over-represented there too.
The reviewer's point stands: to actually escape shared study bias you need an unbiased
substrate. HuRI (Luck et al., Nature 2020) is exactly that -- a proteome-scale, systematic
all-by-all yeast-two-hybrid map, not assembled from the literature, so it does not favour
famous proteins.

The result here is a COVERAGE result, not a signal result, and the distinction is the whole
point: HuRI turns out to be systematically blind to the centrosome module (large coiled-coil
scaffolds express poorly in Y2H). PLK1 itself has HuRI degree 1. So HuRI cannot confirm or
refute Nsp13-PLK1 -- a 0-overlap here means "the assay can't see this biology," not "the
prediction is wrong." Reporting it any other way would repeat the exact error this validation
suite exists to prevent.

Data: HuRI from NDEx (owner ccsb = the Vidal lab), UUID 73bc2c06-5fb2-11e9-9f06-0ac135e8bacf,
parsed once to data/huri_edges.csv (gene-symbol edge list).
"""
import random
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HURI_EDGES = "data/huri_edges.csv"
CENTROSOME_MODULE = ["PLK1", "CEP135", "CEP43", "CNTRL", "NIN", "NINL", "PCNT",
                     "CDK5RAP2", "AKAP9", "CENPF", "CEP250", "CEP68", "CEP112", "CEP350"]
# the 6 Nsp13 preys that ARE IntAct physical partners of PLK1 (test 12's positive evidence)
INTACT_SUPPORT = ["CEP135", "CEP43", "CNTRL", "NIN", "NINL", "TBK1"]
N_NULL = 2000


def load_huri():
    df = pd.read_csv(HURI_EDGES)
    adj = defaultdict(set)
    for a, b in zip(df["geneA"], df["geneB"]):
        adj[a].add(b)
        adj[b].add(a)
    return adj


def main():
    random.seed(20240704)
    adj = load_huri()
    present = set(adj)
    gordon = pd.read_csv("data/interactome.csv")
    print(f"HuRI: {len(present)} genes with >=1 Y2H interaction, {sum(len(v) for v in adj.values())//2} edges\n")

    # 1) The decisive coverage diagnosis
    print("=== Coverage of the centrosome module in HuRI (the module Nsp13 binds) ===")
    cov_rows = []
    for g in CENTROSOME_MODULE:
        d = len(adj.get(g, ()))
        cov_rows.append({"gene": g, "huri_degree": d, "present": g in present})
        flag = "" if g in present else "  <- ABSENT from HuRI"
        print(f"  {g:10s} degree {d}{flag}")
    absent = [r["gene"] for r in cov_rows if not r["present"]]
    print(f"\n  {len(absent)}/{len(CENTROSOME_MODULE)} centrosome-module proteins are entirely absent from HuRI: {absent}")
    print(f"  PLK1 HuRI degree: {len(adj.get('PLK1', ()))} (neighbor(s): {sorted(adj.get('PLK1', ())) or 'none'})")

    # the specific IntAct-supporting evidence -- is any of it even visible to Y2H?
    support_present = [g for g in INTACT_SUPPORT if g in present]
    print(f"\n  Of the 6 Nsp13 preys that support PLK1 in IntAct {INTACT_SUPPORT},")
    print(f"  only {len(support_present)} appear in HuRI at all: {support_present}")

    # 2) The formal overlap test, reported honestly as UNINFORMATIVE
    plk1_partners = adj.get("PLK1", set())
    nsp13_preys = set(gordon.loc[gordon["bait"] == "Nsp13", "prey_gene"])
    real_overlap = nsp13_preys & plk1_partners
    print(f"\n=== Formal overlap test (Nsp13 preys vs PLK1's HuRI partners) ===")
    print(f"  PLK1 has {len(plk1_partners)} HuRI partner(s); Nsp13 prey overlap = {len(real_overlap)} {sorted(real_overlap)}")
    print(f"  VERDICT: UNINFORMATIVE -- with PLK1 at degree {len(plk1_partners)}, no gene set could overlap it. "
          f"This is a coverage failure of Y2H, not evidence about the hypothesis.")

    # 3) Per-bait HuRI coverage, to show the centrosome (Nsp13) is a specific Y2H blind spot
    print(f"\n=== HuRI coverage per bait (fraction of preys present in HuRI, baits with >=5 preys) ===")
    rows = []
    for bait, sub in gordon.groupby("bait"):
        preys = sub["prey_gene"].tolist()
        if len(preys) < 5:
            continue
        pres = [g for g in preys if g in present]
        internal = sum(1 for g in pres for h in adj[g] if h in set(pres)) // 2
        rows.append({"bait": bait, "n_preys": len(preys), "n_in_huri": len(pres),
                     "coverage": round(len(pres) / len(preys), 2), "internal_huri_edges": internal})
    cov = pd.DataFrame(rows).sort_values("coverage")
    print(cov.to_string(index=False))

    cov.to_csv("data/validation/huri_coverage_per_bait.csv", index=False)
    pd.DataFrame(cov_rows).to_csv("data/validation/huri_centrosome_coverage.csv", index=False)

    # figure: HuRI degree of the centrosome module -- most bars at zero => Y2H blind spot
    fig, ax = plt.subplots(figsize=(10, 5.2))
    genes = [r["gene"] for r in cov_rows]
    degs = [r["huri_degree"] for r in cov_rows]
    colors = ["#c0392b" if g == "PLK1" else ("#b9c6d6" if d == 0 else "#3b6fa0")
              for g, d in zip(genes, degs)]
    ax.bar(range(len(genes)), degs, color=colors, edgecolor="white")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_ylabel("HuRI (systematic Y2H) interaction degree")
    ax.set_title("HuRI is blind to the centrosome module Nsp13 binds\n"
                 f"{len(absent)} of {len(genes)} proteins have zero Y2H interactions; PLK1 (red) has 1 "
                 "— so HuRI cannot test Nsp13→PLK1")
    for i, (g, d) in enumerate(zip(genes, degs)):
        if d == 0:
            ax.text(i, 0.3, "absent", rotation=90, ha="center", va="bottom", fontsize=7, color="#7a8899")
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#c0392b", label="PLK1 (the predicted candidate) — degree 1"),
        Patch(facecolor="#3b6fa0", label="present in HuRI"),
        Patch(facecolor="#b9c6d6", label="absent from HuRI (Y2H blind spot)"),
    ], fontsize=8, frameon=False)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig("output/validation/huri_centrosome_blindspot.png", dpi=150)
    print("  saved output/validation/huri_centrosome_blindspot.png")

    print("\n=== Conclusion ===")
    print("  HuRI cannot adjudicate Nsp13-PLK1: PLK1 and most of the centrosome module are")
    print("  effectively invisible to systematic Y2H. Every dataset that DOES cover the centrosome")
    print("  (STRING, IntAct, the Gordon AP-MS itself) is literature-curated -- so the shared-study-")
    print("  bias caveat on the PLK1 result cannot be removed with existing unbiased binary data.")
    print("  Closing that gap requires a targeted wet-lab assay (co-IP / pulldown of Nsp13 + PLK1),")
    print("  not another database. That is the honest ceiling of what computation can establish here.")
    print("\n  saved data/validation/huri_coverage_per_bait.csv, huri_centrosome_coverage.csv")


if __name__ == "__main__":
    main()
