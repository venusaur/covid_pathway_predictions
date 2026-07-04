"""Test 8: external validation via IntAct (EBI) -- a database our pipeline never
touched (not STRING, not the Gordon et al. NDEx network, not Enrichr). IntAct
aggregates curated interactions from many independent SARS-CoV-2 interactome
papers (Gordon, Stukalov, Li, Gao, Chen, Brandherm, Thorne, etc.), each tagged
with its own PubMed ID, so we can check whether an INDEPENDENT dataset
(different paper, different method) supports the specific predicted edge --
not just "this protein binds some SARS-CoV-2 protein somewhere."

Query PLK1 and TOMM70 by UniProt accession (a properly-scoped, reliable query --
a naive free-text search for "nsp13"/"nsp4" against IntAct returned ~984K hits,
i.e. was not correctly scoped, and is NOT used here).
"""
import time
import urllib.request
from io import StringIO
import pandas as pd

PSICQUIC = "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/{}?format=tab25"
COLUMNS = ["idA", "idB", "altA", "altB", "aliasA", "aliasB", "method", "author",
           "pubid", "taxidA", "taxidB", "interactionType", "sourceDb", "interactionId", "confidence"]

SARS_COV_2_TAXID = "2697049"
GORDON_PMID = "32353859"

TARGETS = {
    "TOMM70": {"query": "tomm70", "predicted_bait": "Nsp4", "our_accession": "O94826"},
    "PLK1": {"query": "plk1", "predicted_bait": "Nsp13", "our_accession": "P53350"},
}


def fetch(query):
    url = PSICQUIC.format(query)
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()
    if not text.strip():
        return pd.DataFrame(columns=COLUMNS)
    return pd.read_csv(StringIO(text), sep="\t", header=None, names=COLUMNS)


def uniprot_id(field):
    return field.split(":")[1] if ":" in field else field


def main():
    all_rows = []
    for gene, cfg in TARGETS.items():
        print(f"\n=== {gene} (predicted interactor of {cfg['predicted_bait']}) ===")
        df = fetch(cfg["query"])
        print(f"  total IntAct records: {len(df)}")

        viral = df[df["taxidA"].str.contains(SARS_COV_2_TAXID) | df["taxidB"].str.contains(SARS_COV_2_TAXID)]
        print(f"  records involving a SARS-CoV-2 partner: {len(viral)}")

        # identify which UniProt accession the viral partner actually is
        viral_partners = set()
        for _, row in viral.iterrows():
            a_is_viral = SARS_COV_2_TAXID in row["taxidA"]
            partner = uniprot_id(row["idA"] if a_is_viral else row["idB"])
            viral_partners.add(partner)
        print(f"  distinct SARS-CoV-2 UniProt partner(s): {viral_partners}")

        pmids = set()
        for p in viral["pubid"]:
            for token in str(p).split("|"):
                if token.startswith("pubmed:"):
                    pmids.add(token.replace("pubmed:", ""))
        independent_pmids = pmids - {GORDON_PMID}
        print(f"  publications reporting this (excl. Gordon et al. {GORDON_PMID}): {sorted(independent_pmids)}")

        all_rows.append({
            "candidate_gene": gene,
            "predicted_bait": cfg["predicted_bait"],
            "total_intact_records": len(df),
            "sars_cov2_records": len(viral),
            "sars_cov2_uniprot_partners": ";".join(sorted(viral_partners)),
            "independent_pmids_confirming_ANY_sars_cov2_partner": ";".join(sorted(independent_pmids)),
            "confirms_predicted_bait_specifically": cfg["predicted_bait"] in str(viral_partners),
        })
        time.sleep(1)

    out = pd.DataFrame(all_rows)
    out.to_csv("data/validation/external_validation_intact.csv", index=False)
    print("\n=== SUMMARY ===")
    print(out.to_string(index=False))
    print("\nsaved data/validation/external_validation_intact.csv")


if __name__ == "__main__":
    main()
