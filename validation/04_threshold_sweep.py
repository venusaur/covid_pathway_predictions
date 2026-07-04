"""Test 6: sweep the researcher-chosen parameters behind the PLK1/TOMM70 calls
(STRING required_score confidence threshold, and add_nodes expansion width).
A real result should be stable across a reasonable band; a p-hacked one lives
at one specific setting.
"""
import time
import urllib.request
import urllib.parse
from io import StringIO
import pandas as pd

STRING_URL = "https://string-db.org/api/tsv/network"

_df = pd.read_csv("data/interactome.csv")


def _prey_list(bait):
    # pulled programmatically from the real data, not hand-typed, to avoid transcription errors
    return sorted(_df.loc[_df["bait"] == bait, "prey_gene"].tolist())


CASES = {
    "Nsp13": {"genes": _prey_list("Nsp13"), "target": "PLK1"},
    "Nsp4": {"genes": _prey_list("Nsp4"), "target": "TOMM70"},
}
for bait, cfg in CASES.items():
    print(f"{bait} full prey list (n={len(cfg['genes'])}): {cfg['genes']}")

REQUIRED_SCORES = [400, 500, 600, 700, 800, 900]  # STRING scale 0-1000 = 0.0-1.0 confidence
ADD_NODES_VALUES = [5, 10, 15, 20]


def string_network(genes, add_nodes=10, required_score=None):
    ids = "%0d".join(genes)
    params = {"identifiers": ids, "species": "9606", "add_nodes": add_nodes}
    if required_score is not None:
        params["required_score"] = required_score
    url = STRING_URL + "?" + urllib.parse.urlencode(params, safe="%0d")
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()
    return pd.read_csv(StringIO(text), sep="\t")


def added_nodes(edges, known):
    return (set(edges["preferredName_A"]) | set(edges["preferredName_B"])) - known


def main():
    rows = []
    for bait, cfg in CASES.items():
        known = set(cfg["genes"])
        target = cfg["target"]

        print(f"\n=== {bait}: sweeping required_score (add_nodes=10) ===")
        for rs in REQUIRED_SCORES:
            try:
                edges = string_network(cfg["genes"], add_nodes=10, required_score=rs)
                found = target in added_nodes(edges, known)
            except Exception as ex:
                found = f"ERROR: {ex}"
            print(f"  required_score={rs/1000:.1f}: {target} predicted = {found}")
            rows.append({"bait": bait, "target": target, "sweep": "required_score",
                         "value": rs/1000, "predicted": found})
            time.sleep(1)

        print(f"=== {bait}: sweeping add_nodes width (default score threshold) ===")
        for an in ADD_NODES_VALUES:
            try:
                edges = string_network(cfg["genes"], add_nodes=an)
                found = target in added_nodes(edges, known)
            except Exception as ex:
                found = f"ERROR: {ex}"
            print(f"  add_nodes={an}: {target} predicted = {found}")
            rows.append({"bait": bait, "target": target, "sweep": "add_nodes",
                         "value": an, "predicted": found})
            time.sleep(1)

    out = pd.DataFrame(rows)
    out.to_csv("data/validation/threshold_sweep.csv", index=False)
    print("\nsaved data/validation/threshold_sweep.csv")

    print("\n=== Summary: fraction of settings where target was predicted ===")
    summary = out.groupby(["bait", "sweep"])["predicted"].apply(lambda s: (s == True).mean())
    print(summary.to_string())


if __name__ == "__main__":
    main()
