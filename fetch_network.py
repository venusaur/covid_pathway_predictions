"""Fetch and parse the real Gordon et al. 2020 SARS-CoV-2-human interactome from NDEx.

Source: kroganlab HEK293T_SARS-CoV-2 network, NDEx UUID 43803262-6d69-11ea-bfdc-0ac135e8bacf
(the 332 high-confidence PPIs behind Fig. 3 of the paper).
"""
import json
import urllib.request
import pandas as pd

NDEX_UUID = "43803262-6d69-11ea-bfdc-0ac135e8bacf"
NDEX_URL = f"https://www.ndexbio.org/v2/network/{NDEX_UUID}"


def fetch_cx():
    with urllib.request.urlopen(NDEX_URL) as r:
        return json.load(r)


def parse(cx):
    aspects = {}
    for asp in cx:
        key = list(asp.keys())[0]
        aspects.setdefault(key, []).extend(asp[key])

    node_name = {n["@id"]: n["n"] for n in aspects["nodes"]}
    is_bait = {}
    for a in aspects.get("nodeAttributes", []):
        if a["n"] == "Bait_Boolean":
            is_bait[a["po"]] = a["v"] == "1"

    edge_attrs = {}
    for a in aspects.get("edgeAttributes", []):
        edge_attrs.setdefault(a["po"], {})[a["n"]] = a["v"]

    rows = []
    for e in aspects["edges"]:
        eid, s, t = e["@id"], e["s"], e["t"]
        attrs = edge_attrs.get(eid, {})
        s_bait, t_bait = is_bait.get(s, False), is_bait.get(t, False)
        # keep only viral(bait)-human(prey) interactions, skip curated human-human edges
        if s_bait == t_bait:
            continue
        bait_id, prey_id = (s, t) if s_bait else (t, s)
        rows.append({
            "bait": node_name[bait_id],
            "prey_gene": node_name[prey_id],
            "prey_uniprot": attrs.get("Preys"),
            "MIST": float(attrs["MIST"]) if attrs.get("MIST") not in (None, "NaN") else None,
            "BFDR": float(attrs["BFDR"]) if attrs.get("BFDR") not in (None, "NaN") else None,
            "AvgSpec": float(attrs["AvgSpec"]) if attrs.get("AvgSpec") not in (None, "NaN") else None,
            "FoldChange": attrs.get("FoldChange"),
            "Complex": attrs.get("Complex"),
        })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = parse(fetch_cx())
    df.to_csv("data/interactome.csv", index=False)
    print(f"{len(df)} viral-host interactions across {df['bait'].nunique()} baits")
    print(df["bait"].value_counts())
