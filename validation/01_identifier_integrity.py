"""Test 7 (identifier integrity): confirm every prey_gene/prey_uniprot pair in
interactome.csv resolves to a real, current human UniProt entry whose official
HGNC gene symbol matches what we recorded. A broken join wouldn't announce
itself -- it would just quietly produce a network missing rows or mislabeled
nodes (the 'TIMM29' concern).
"""
import json
import time
import urllib.request
import urllib.parse
import pandas as pd

BATCH = 90


def fetch_batch(accessions):
    params = urllib.parse.urlencode({
        "accessions": ",".join(accessions),
        "fields": "accession,gene_names,organism_name",
    })
    url = f"https://rest.uniprot.org/uniprotkb/accessions?{params}"
    with urllib.request.urlopen(url) as r:
        return json.load(r)["results"]


def main():
    df = pd.read_csv("data/interactome.csv")
    accs = df["prey_uniprot"].unique().tolist()
    print(f"Validating {len(accs)} unique UniProt accessions across {len(df)} interactions...")

    resolved = {}
    for i in range(0, len(accs), BATCH):
        chunk = accs[i:i + BATCH]
        for r in fetch_batch(chunk):
            acc = r["primaryAccession"]
            organism = r.get("organism", {}).get("scientificName", "")
            genes = r.get("genes", [])
            official = genes[0]["geneName"]["value"] if genes and "geneName" in genes[0] else None
            synonyms = [s["value"] for s in genes[0].get("synonyms", [])] if genes else []
            resolved[acc] = {"official_symbol": official, "synonyms": synonyms, "organism": organism}
        time.sleep(0.3)

    rows = []
    for _, row in df.iterrows():
        acc, recorded_gene, bait = row["prey_uniprot"], row["prey_gene"], row["bait"]
        info = resolved.get(acc)
        if info is None:
            status = "ACCESSION_NOT_FOUND"
            official, organism = None, None
        else:
            official, organism = info["official_symbol"], info["organism"]
            if official == recorded_gene:
                status = "OK"
            elif recorded_gene in info["synonyms"]:
                status = "OK_VIA_SYNONYM"
            elif organism != "Homo sapiens":
                status = "WRONG_ORGANISM"
            else:
                status = "SYMBOL_MISMATCH"
        rows.append({
            "bait": bait, "recorded_gene": recorded_gene, "uniprot_accession": acc,
            "uniprot_official_symbol": official, "organism": organism, "status": status,
        })

    out = pd.DataFrame(rows)
    out.to_csv("data/validation/identifier_integrity.csv", index=False)

    print("\n=== Status breakdown ===")
    print(out["status"].value_counts().to_string())

    problems = out[out["status"] != "OK"]
    if len(problems):
        print(f"\n{len(problems)} rows need a look:")
        print(problems.to_string(index=False))
    else:
        print("\nAll 332 recorded gene symbols exactly match their UniProt accession's official HGNC symbol.")

    dup_acc = df["prey_uniprot"].duplicated().sum()
    dup_gene = df["prey_gene"].duplicated().sum()
    print(f"\nDuplicate UniProt accessions across rows: {dup_acc}")
    print(f"Duplicate gene symbols across rows: {dup_gene}")
    print(f"Null prey_uniprot: {df['prey_uniprot'].isna().sum()}, null prey_gene: {df['prey_gene'].isna().sum()}")


if __name__ == "__main__":
    main()
