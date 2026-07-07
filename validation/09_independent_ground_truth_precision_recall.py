"""Next-step 2: re-run precision/recall against a GENUINELY independent ground truth.

Test 7 (06_holdout_precision_recall.py) only ever split Gordon et al.'s own 332
interactions into train/test folds -- a real signal-vs-null comparison, but still
internal to one dataset. This script is a stronger test: take Gordon et al.'s FULL real
prey list for a bait (no holdout needed) as input to the identical STRING
guilt-by-association procedure, and check how many candidates match interactors that
Stukalov et al. 2021 -- a separate lab, separate cell line, separate AP-MS pipeline --
independently reported for the SAME viral protein but that Gordon's own screen missed.

Recall here is defined against "novel Stukalov interactors": Stukalov's reported hits
for a bait, minus whatever Gordon et al. already found (recovering those would be
trivial, since they're literally the query input). Compared against a degree-matched
null of the same construction as validation/05_null_model.py and validation/06.
"""
import time
import random
import urllib.request
import urllib.parse
from io import StringIO
import openpyxl
import pandas as pd

STRING_URL = "https://string-db.org/api/tsv/network"
ADD_NODES = 15
MIN_GORDON_PREYS = 3
STUKALOV_XLSX = "data/validation/stukalov_supp_table2.xlsx"

# Gordon et al. bait name -> Stukalov et al. bait name, for the baits both studies tested
BAIT_MAP = {
    "E": "E", "M": "M", "N": "N", "Nsp1": "NSP1", "Nsp12": "NSP12", "Nsp13": "NSP13",
    "Nsp14": "NSP14", "Nsp15": "NSP15", "Nsp2": "NSP2", "Nsp4": "NSP4", "Nsp6": "NSP6",
    "Nsp7": "NSP7", "Nsp8": "NSP8", "Nsp9": "NSP9", "Orf3a": "ORF3", "Orf6": "ORF6",
    "Orf7a": "ORF7a", "Orf8": "ORF8", "Orf9b": "ORF9b", "Spike": "S",
}


def string_network(genes, add_nodes=ADD_NODES):
    ids = "%0d".join(genes)
    params = {"identifiers": ids, "species": "9606", "add_nodes": add_nodes}
    url = STRING_URL + "?" + urllib.parse.urlencode(params, safe="%0d")
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()
    return pd.read_csv(StringIO(text), sep="\t")


def load_stukalov():
    wb = openpyxl.load_workbook(STUKALOV_XLSX, read_only=True)
    ws = wb["A - Significant interactions"]
    rows = ws.iter_rows(values_only=True)
    header = next(rows)
    df = pd.DataFrame(rows, columns=header)
    return df[df["bait_organism"] == "SARS-CoV-2"]


def one_trial(query_genes, target_genes):
    if not target_genes:
        return None
    try:
        edges = string_network(query_genes)
        candidates = (set(edges["preferredName_A"]) | set(edges["preferredName_B"])) - set(query_genes)
    except Exception as ex:
        print(f"    [error] {ex}")
        return None
    hit = candidates & target_genes
    recall = len(hit) / len(target_genes)
    precision = len(hit) / len(candidates) if candidates else float("nan")
    return {"n_candidates": len(candidates), "n_targets": len(target_genes),
            "n_recovered": len(hit), "recall": recall, "precision": precision,
            "recovered_genes": ";".join(sorted(hit))}


def main():
    random.seed(20240704)
    gordon = pd.read_csv("data/interactome.csv")
    pool = gordon["prey_gene"].tolist()
    stukalov = load_stukalov()

    rows = []
    for gordon_bait, stukalov_bait in BAIT_MAP.items():
        gordon_genes = gordon.loc[gordon["bait"] == gordon_bait, "prey_gene"].tolist()
        if len(gordon_genes) < MIN_GORDON_PREYS:
            print(f"{gordon_bait}: skipped, only {len(gordon_genes)} Gordon preys (< {MIN_GORDON_PREYS})")
            continue

        stukalov_genes = set(stukalov.loc[stukalov["bait_name"] == stukalov_bait, "gene_name"].unique())
        novel_targets = stukalov_genes - set(gordon_genes)

        print(f"\n=== {gordon_bait} (Gordon n={len(gordon_genes)}) vs Stukalov's {stukalov_bait} "
              f"(n={len(stukalov_genes)}, {len(novel_targets)} novel beyond Gordon's own list) ===")
        if not novel_targets:
            print("  no novel Stukalov targets to test recall against -- skipping")
            continue

        real = one_trial(gordon_genes, novel_targets)
        if real:
            real.update({"bait": gordon_bait, "condition": "real"})
            rows.append(real)
            print(f"  [real] recall={real['recall']:.2f} precision={real['precision']:.2f} "
                  f"({real['n_recovered']}/{real['n_targets']} recovered)"
                  + (f" -- hit: {real['recovered_genes']}" if real['n_recovered'] else ""))
        time.sleep(1)

        null_genes = random.sample(pool, len(gordon_genes))
        null_res = one_trial(null_genes, novel_targets)
        if null_res:
            null_res.update({"bait": gordon_bait, "condition": "null_baseline"})
            rows.append(null_res)
            print(f"  [null] recall={null_res['recall']:.2f} precision={null_res['precision']:.2f} "
                  f"({null_res['n_recovered']}/{null_res['n_targets']} recovered)")
        time.sleep(1)

    out = pd.DataFrame(rows)
    out.to_csv("data/validation/independent_ground_truth_precision_recall.csv", index=False)

    print("\n=== SUMMARY: mean recall / precision, real (vs Stukalov ground truth) vs degree-matched null ===")
    print(out.groupby("condition")[["recall", "precision"]].mean().to_string())

    print("\n=== Per-bait recall, real vs null ===")
    print(out.groupby(["bait", "condition"])["recall"].mean().unstack().to_string())

    print(f"\nBaits actually testable (had >=1 novel Stukalov target): "
          f"{sorted(out['bait'].unique())}")
    print("saved data/validation/independent_ground_truth_precision_recall.csv")


if __name__ == "__main__":
    main()
