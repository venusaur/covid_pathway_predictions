"""Next-step 5: does Nsp13->PLK1 (and Nsp4->TOMM70) reproduce under a STRUCTURALLY
DIFFERENT network substrate?

The original predictions came from STRING, whose combined score blends seven evidence
channels -- and test 11 showed the method's successes are driven largely by STRING's
text-mining channel (literature co-mention). A genuinely independent check is to redo the
guilt-by-association on a network built a completely different way: IntAct's curated
*physical* interactions (experimentally-observed binding from primary literature, no
co-expression or text-mining inference). If PLK1 still emerges for Nsp13 on IntAct's
graph, that's cross-method reproduction; if it evaporates, the call was STRING-specific.

Two tests:
  A) Overlap significance -- are the bait's real preys enriched among the candidate's
     IntAct physical partners, versus a degree-matched null drawn from the 332-gene pool?
  B) Substrate-swapped re-ranking -- redo the missing-component ranking using IntAct
     partners instead of STRING add_nodes, and report where the candidate lands.
"""
import time
import re
import random
import urllib.request
from io import StringIO
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PSICQUIC = ("https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/"
            "webservices/current/search/query/{}?format=tab25")
COLS = ["idA", "idB", "altA", "altB", "aliasA", "aliasB", "method", "author",
        "pubid", "taxidA", "taxidB", "itype", "src", "iid", "conf"]
N_NULL = 2000
GENE_RE = re.compile(r"uniprotkb:([A-Za-z0-9_.-]+)\(gene name\)")

_cache = {}


def _gene(alias):
    m = GENE_RE.search(str(alias))
    return m.group(1) if m else None


def intact_human_partners(accession, self_gene):
    """Distinct human gene symbols that physically interact with `accession` in IntAct."""
    if accession in _cache:
        return _cache[accession]
    url = PSICQUIC.format(accession)
    try:
        text = urllib.request.urlopen(url).read().decode()
    except Exception as ex:
        print(f"    [error {accession}] {ex}")
        _cache[accession] = set()
        return set()
    if not text.strip():
        _cache[accession] = set()
        return set()
    df = pd.read_csv(StringIO(text), sep="\t", header=None, names=COLS)
    hh = df[df["taxidA"].str.contains("9606", na=False) & df["taxidB"].str.contains("9606", na=False)]
    partners = set()
    for _, r in hh.iterrows():
        for g in (_gene(r["aliasA"]), _gene(r["aliasB"])):
            if g and g != self_gene:
                partners.add(g)
    _cache[accession] = partners
    time.sleep(1)
    return partners


def overlap_test(bait, bait_preys, candidate_gene, candidate_acc, gordon):
    random.seed(20240704)
    partners = intact_human_partners(candidate_acc, candidate_gene)
    prey_set = set(bait_preys)
    real_overlap = prey_set & partners
    pool = gordon["prey_gene"].tolist()
    null_overlaps = []
    for _ in range(N_NULL):
        s = set(random.sample(pool, len(prey_set)))
        null_overlaps.append(len(s & partners))
    real = len(real_overlap)
    ge = sum(1 for x in null_overlaps if x >= real)
    p = (ge + 1) / (N_NULL + 1)
    mean_null = sum(null_overlaps) / len(null_overlaps)
    return {
        "bait": bait, "candidate": candidate_gene,
        "n_candidate_intact_partners": len(partners),
        "n_bait_preys": len(prey_set),
        "real_overlap": real,
        "overlap_genes": ";".join(sorted(real_overlap)),
        "null_mean_overlap": round(mean_null, 3),
        "empirical_p": round(p, 4),
        "_null_dist": null_overlaps,
    }


def plot_reproduction(results):
    fig, axes = plt.subplots(1, len(results), figsize=(6.4 * len(results), 4.8))
    if len(results) == 1:
        axes = [axes]
    for ax, res in zip(axes, results):
        null = res["_null_dist"]
        real = res["real_overlap"]
        maxx = max(max(null), real) + 1
        bins = range(0, maxx + 1)
        ax.hist(null, bins=bins, align="left", color="#b9c6d6", edgecolor="white",
                label=f"null (random {res['n_bait_preys']}-gene sets)")
        ax.axvline(real, color="#b23a2e", lw=2.5, zorder=5,
                   label=f"real {res['bait']} preys = {real}")
        ax.set_xlabel(f"# of {res['bait']}'s preys that are\nIntAct physical partners of {res['candidate']}")
        ax.set_ylabel("null draws")
        verdict = "reproduces" if res["empirical_p"] < 0.05 else "does not reproduce"
        ax.set_title(f"{res['bait']} → {res['candidate']} on IntAct: {verdict}\n"
                     f"empirical p = {res['empirical_p']}  "
                     f"({res['n_candidate_intact_partners']} {res['candidate']} partners total)",
                     fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        ax.grid(True, axis="y", alpha=0.25)
    fig.suptitle("Cross-method check: does the prediction survive swapping STRING's blended score\n"
                 "for IntAct's curated physical-interaction graph?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("output/validation/alt_substrate_reproduction.png", dpi=150)
    print("saved output/validation/alt_substrate_reproduction.png")


def substrate_swapped_ranking(bait, gordon, target_gene):
    """IntAct analog of STRING add_nodes: rank human proteins by how many of the bait's
    preys they physically partner with in IntAct."""
    sub = gordon[gordon["bait"] == bait]
    prey_genes = set(sub["prey_gene"])
    counts = {}
    for _, row in sub.iterrows():
        partners = intact_human_partners(row["prey_uniprot"], row["prey_gene"])
        for p in partners:
            if p not in prey_genes:  # candidate = not already a known prey
                counts[p] = counts.get(p, 0) + 1
    ranked = sorted(counts.items(), key=lambda x: -x[1])
    rank = next((i + 1 for i, (g, _) in enumerate(ranked) if g == target_gene), None)
    top = ranked[:10]
    return {
        "bait": bait, "target": target_gene,
        "target_rank": rank,
        "target_edge_count": counts.get(target_gene, 0),
        "n_total_candidates": len(ranked),
        "top10": ";".join(f"{g}({c})" for g, c in top),
    }


def main():
    gordon = pd.read_csv("data/interactome.csv")
    # candidate accessions (not preys themselves): PLK1, TOMM70
    CAND = {"PLK1": "P53350", "TOMM70": "O94826"}

    print("=== TEST A: overlap of bait preys with candidate's IntAct physical partners ===")
    rowsA = []
    for bait, cand in [("Nsp13", "PLK1"), ("Nsp4", "TOMM70")]:
        preys = gordon.loc[gordon["bait"] == bait, "prey_gene"].tolist()
        res = overlap_test(bait, preys, cand, CAND[cand], gordon)
        rowsA.append(res)
        print(f"  {bait}->{cand}: {res['real_overlap']}/{res['n_bait_preys']} preys are IntAct "
              f"partners of {cand} (null mean {res['null_mean_overlap']}), "
              f"empirical p={res['empirical_p']}"
              + (f" -- overlap: {res['overlap_genes']}" if res['real_overlap'] else ""))
    plot_reproduction(rowsA)
    pd.DataFrame([{k: v for k, v in r.items() if k != "_null_dist"} for r in rowsA]).to_csv(
        "data/validation/alt_substrate_overlap.csv", index=False)

    print("\n=== TEST B: substrate-swapped re-ranking (IntAct instead of STRING) ===")
    rowsB = []
    for bait, target in [("Nsp13", "PLK1"), ("Nsp4", "TOMM70")]:
        res = substrate_swapped_ranking(bait, gordon, target)
        rowsB.append(res)
        rank_str = f"#{res['target_rank']}" if res['target_rank'] else "NOT RANKED"
        print(f"  {bait}: {target} ranks {rank_str} of {res['n_total_candidates']} IntAct candidates "
              f"(connects to {res['target_edge_count']} of the bait's preys)")
        print(f"     top IntAct candidates: {res['top10']}")
    pd.DataFrame(rowsB).to_csv("data/validation/alt_substrate_ranking.csv", index=False)

    print("\nsaved data/validation/alt_substrate_overlap.csv and alt_substrate_ranking.csv")


if __name__ == "__main__":
    main()
