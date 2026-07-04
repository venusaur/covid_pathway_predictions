"""Test 1: degree-preserving null model.

Key structural fact checked first: in the real 332-interaction dataset, every
prey gene is unique to exactly one bait (prey-side degree is always 1 -- see
identifier_integrity.csv / no duplicate prey_gene across the whole table).
That means a bipartite double-edge-swap configuration model, which must
preserve each node's degree exactly, degenerates here to: randomly partition
the fixed pool of 332 prey genes into groups of the SAME SIZES as the real
bait groups. This is provably degree-preserving on both sides for this
dataset, so we use it directly rather than running a swap chain that would
converge to the same distribution.

For each shuffle we draw a random gene set the same size as Nsp13's (40) and
Nsp4's (8) real prey count, feed it through the identical STRING add_nodes
prediction, and check how often PLK1 / TOMM70 (or any of STRING's best-known
hub proteins) surface under labels that carry no real biology. If they surface
about as often under randomized labels as under the true labels, the original
call is a database/hub artifact, not signal.
"""
import time
import urllib.request
import urllib.parse
from io import StringIO
import random
import pandas as pd

STRING_URL = "https://string-db.org/api/tsv/network"
N_SHUFFLES = 25


def string_network(genes, add_nodes=10):
    ids = "%0d".join(genes)
    params = {"identifiers": ids, "species": "9606", "add_nodes": add_nodes}
    url = STRING_URL + "?" + urllib.parse.urlencode(params, safe="%0d")
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()
    return pd.read_csv(StringIO(text), sep="\t")


def run_null(pool, group_size, target, label, n=N_SHUFFLES):
    hits = 0
    top_candidates = []
    for i in range(n):
        genes = random.sample(pool, group_size)
        try:
            edges = string_network(genes, add_nodes=10)
            added = (set(edges["preferredName_A"]) | set(edges["preferredName_B"])) - set(genes)
            hit = target in added
            top_candidates.extend(added)
        except Exception as ex:
            hit = False
            print(f"  [error on shuffle {i}] {ex}")
        hits += hit
        print(f"  shuffle {i+1}/{n}: {target} present = {hit}")
        time.sleep(1)
    return hits, top_candidates


def main():
    random.seed(20240704)
    df = pd.read_csv("data/interactome.csv")
    pool = df["prey_gene"].tolist()  # 332 unique genes, each degree 1 -- see docstring
    assert df["prey_gene"].duplicated().sum() == 0, "assumption violated: a prey appears under >1 bait"

    print(f"Pool size: {len(pool)} (each gene degree-1 by construction)")

    print(f"\n=== Null for 'Nsp13'-sized (n=40) random gene sets: does PLK1 surface? ===")
    plk1_hits, plk1_candidates = run_null(pool, 40, "PLK1", "Nsp13-null")

    print(f"\n=== Null for 'Nsp4'-sized (n=8) random gene sets: does TOMM70 surface? ===")
    tomm70_hits, tomm70_candidates = run_null(pool, 8, "TOMM70", "Nsp4-null")

    from collections import Counter
    plk1_rate = plk1_hits / N_SHUFFLES
    tomm70_rate = tomm70_hits / N_SHUFFLES

    print(f"\n=== RESULTS ===")
    print(f"PLK1 predicted under {N_SHUFFLES} random 40-gene null sets: {plk1_hits}/{N_SHUFFLES} ({plk1_rate:.1%})")
    print(f"  (real Nsp13 labels: PLK1 predicted = True, per validation/04_threshold_sweep.py)")
    print(f"TOMM70 predicted under {N_SHUFFLES} random 8-gene null sets: {tomm70_hits}/{N_SHUFFLES} ({tomm70_rate:.1%})")
    print(f"  (real Nsp4 labels: TOMM70 predicted = True only at loosest threshold, per sweep)")

    print("\nMost common candidates surfacing under randomized 40-gene null (top 15):")
    for gene, count in Counter(plk1_candidates).most_common(15):
        print(f"  {gene}: {count}/{N_SHUFFLES}")

    print("\nMost common candidates surfacing under randomized 8-gene null (top 15):")
    for gene, count in Counter(tomm70_candidates).most_common(15):
        print(f"  {gene}: {count}/{N_SHUFFLES}")

    out = pd.DataFrame([
        {"bait": "Nsp13", "target": "PLK1", "group_size": 40, "n_shuffles": N_SHUFFLES,
         "null_hit_rate": plk1_rate, "real_label_predicted": True},
        {"bait": "Nsp4", "target": "TOMM70", "group_size": 8, "n_shuffles": N_SHUFFLES,
         "null_hit_rate": tomm70_rate, "real_label_predicted": True},
    ])
    out.to_csv("data/validation/null_model_summary.csv", index=False)

    pd.DataFrame(Counter(plk1_candidates).most_common(), columns=["gene", "count"]).to_csv(
        "data/validation/null_model_nsp13_candidate_freq.csv", index=False)
    pd.DataFrame(Counter(tomm70_candidates).most_common(), columns=["gene", "count"]).to_csv(
        "data/validation/null_model_nsp4_candidate_freq.csv", index=False)
    print("\nsaved data/validation/null_model_*.csv")


if __name__ == "__main__":
    main()
