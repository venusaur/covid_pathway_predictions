"""Test 3: precision/recall via hold-out, benchmarked against the degree-matched null.

For each bait with >=10 real preys: hide ~20% of its known preys, run the exact
guilt-by-association procedure (STRING add_nodes) on the remaining 80%, and
check how many of the hidden true preys get recovered among the predicted
candidates (recall), and what fraction of predicted candidates are hidden true
positives (precision -- a conservative lower bound, since a "wrong" candidate
could still be a real undiscovered interactor, not necessarily junk).

Repeated over several random holdout splits per bait, and compared against an
identical procedure run on random gene sets of matched size drawn from the
same 332-gene pool (the same null construction as validation/05_null_model.py)
to see whether the real biology actually helps over the degree-matched null.
"""
import time
import random
import urllib.request
import urllib.parse
from io import StringIO
import pandas as pd

STRING_URL = "https://string-db.org/api/tsv/network"
N_ITER = 5
HOLDOUT_FRAC = 0.2
ADD_NODES = 15
MIN_PREYS = 10


def string_network(genes, add_nodes=ADD_NODES):
    ids = "%0d".join(genes)
    params = {"identifiers": ids, "species": "9606", "add_nodes": add_nodes}
    url = STRING_URL + "?" + urllib.parse.urlencode(params, safe="%0d")
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()
    return pd.read_csv(StringIO(text), sep="\t")


def one_trial(all_genes, held_out, training):
    try:
        edges = string_network(training)
        candidates = (set(edges["preferredName_A"]) | set(edges["preferredName_B"])) - set(training)
    except Exception as ex:
        print(f"    [error] {ex}")
        return None
    hit = candidates & set(held_out)
    recall = len(hit) / len(held_out) if held_out else float("nan")
    precision = len(hit) / len(candidates) if candidates else float("nan")
    return {"n_candidates": len(candidates), "n_held_out": len(held_out),
            "n_recovered": len(hit), "recall": recall, "precision": precision}


def main():
    random.seed(20240704)
    df = pd.read_csv("data/interactome.csv")
    pool = df["prey_gene"].tolist()
    bait_sizes = df["bait"].value_counts()
    eligible = bait_sizes[bait_sizes >= MIN_PREYS].index.tolist()
    print(f"Eligible baits (>= {MIN_PREYS} preys): {eligible}")

    rows = []
    for bait in eligible:
        genes = df.loc[df["bait"] == bait, "prey_gene"].tolist()
        n_hold = max(1, int(round(len(genes) * HOLDOUT_FRAC)))
        print(f"\n=== {bait} (n={len(genes)}, holding out {n_hold} per iteration) ===")

        for it in range(N_ITER):
            shuffled = genes[:]
            random.shuffle(shuffled)
            held_out, training = shuffled[:n_hold], shuffled[n_hold:]
            res = one_trial(genes, held_out, training)
            if res:
                res.update({"bait": bait, "iter": it, "condition": "real"})
                rows.append(res)
                print(f"  [real] iter {it}: recall={res['recall']:.2f} precision={res['precision']:.2f} "
                      f"({res['n_recovered']}/{res['n_held_out']} recovered, {res['n_candidates']} candidates)")
            time.sleep(1)

            # degree-matched null: same-size random draw from the full pool, same holdout mechanics
            null_group = random.sample(pool, len(genes))
            null_held_out, null_training = null_group[:n_hold], null_group[n_hold:]
            null_res = one_trial(null_group, null_held_out, null_training)
            if null_res:
                null_res.update({"bait": bait, "iter": it, "condition": "null"})
                rows.append(null_res)
                print(f"  [null] iter {it}: recall={null_res['recall']:.2f} precision={null_res['precision']:.2f} "
                      f"({null_res['n_recovered']}/{null_res['n_held_out']} recovered, {null_res['n_candidates']} candidates)")
            time.sleep(1)

    out = pd.DataFrame(rows)
    out.to_csv("data/validation/holdout_precision_recall.csv", index=False)

    print("\n=== SUMMARY: mean recall / precision, real vs degree-matched null ===")
    summary = out.groupby("condition")[["recall", "precision"]].mean()
    print(summary.to_string())

    print("\n=== Per-bait mean recall, real vs null ===")
    per_bait = out.groupby(["bait", "condition"])["recall"].mean().unstack()
    print(per_bait.to_string())

    out.to_csv("data/validation/holdout_precision_recall.csv", index=False)
    print("\nsaved data/validation/holdout_precision_recall.csv")


if __name__ == "__main__":
    main()
