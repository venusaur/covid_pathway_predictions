"""Predict 'missing' host-interactome components per SARS-CoV-2 viral protein.

For each bait's known prey set, ask STRING to expand the network with the
N most-connected additional proteins (guilt-by-association). Any added protein
not already in our 332-interaction dataset is a candidate that plausibly also
binds that viral protein, ranked by STRING confidence to the known preys.
"""
import time
import urllib.request
import urllib.parse
from io import StringIO
import pandas as pd

MIN_PREYS = 5
ADD_NODES = 10
STRING_URL = "https://string-db.org/api/tsv/network"


def string_network(genes, add_nodes=ADD_NODES):
    ids = "%0d".join(genes)
    params = {"identifiers": ids, "species": "9606", "add_nodes": add_nodes}
    url = STRING_URL + "?" + urllib.parse.urlencode(params, safe="%0d")
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()
    return pd.read_csv(StringIO(text), sep="\t")


def main():
    df = pd.read_csv("data/interactome.csv")
    all_preys_by_bait = df.groupby("bait")["prey_gene"].apply(lambda s: set(s.dropna()))

    rows = []
    for bait, genes in all_preys_by_bait.items():
        genes = sorted(genes)
        if len(genes) < MIN_PREYS:
            continue
        print(f"Expanding network for {bait} ({len(genes)} known preys)...")
        try:
            edges = string_network(genes)
        except Exception as ex:
            print(f"  [skip] {ex}")
            continue

        known = set(genes)
        candidates = {}
        for _, e in edges.iterrows():
            a, b, score = e["preferredName_A"], e["preferredName_B"], e["score"]
            for cand, other in [(a, b), (b, a)]:
                if cand not in known and other in known:
                    d = candidates.setdefault(cand, {"scores": [], "linked_to": set()})
                    d["scores"].append(score)
                    d["linked_to"].add(other)

        for cand, d in candidates.items():
            other_baits = sorted(set(df.loc[df["prey_gene"] == cand, "bait"]) - {bait})
            rows.append({
                "bait": bait,
                "candidate_gene": cand,
                "n_edges_to_known_preys": len(d["linked_to"]),
                "max_string_score": max(d["scores"]),
                "connected_to": ";".join(sorted(d["linked_to"])),
                "already_interacts_with_other_baits": ";".join(other_baits) if other_baits else "",
            })
        time.sleep(1)  # be polite to STRING's API

    out = pd.DataFrame(rows).sort_values(
        ["bait", "n_edges_to_known_preys", "max_string_score"], ascending=[True, False, False]
    )
    out.to_csv("data/predicted_missing_components.csv", index=False)
    print(f"\nSaved {len(out)} candidate missing-component predictions across "
          f"{out['bait'].nunique()} baits -> data/predicted_missing_components.csv")


if __name__ == "__main__":
    main()
