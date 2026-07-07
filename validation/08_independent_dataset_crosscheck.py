"""Next-step follow-up: does an INDEPENDENT SARS-CoV-2 AP-MS dataset (not Gordon et al.,
not the data our whole pipeline is built on) list PLK1 among Nsp13's real interactors?

Test 8 (validation/07_external_validation.py) checked IntAct in general and found no
PLK1-Nsp13 record. This script goes one step further and pulls the actual, authoritative
source directly: Stukalov et al. 2021, "Multilevel proteomics reveals host perturbations
by SARS-CoV-2 and SARS-CoV" (Nature 594:246-252, PMID 33845483), Supplementary Table 2
("A - Significant interactions" sheet) -- their own full, curated AP-MS hit list, fetched
straight from Nature's static content server rather than through a third-party database's
partial curation.
"""
import urllib.request
import openpyxl
import pandas as pd

SUPP_TABLE_URL = ("https://static-content.springer.com/esm/"
                   "art%3A10.1038%2Fs41586-021-03493-4/MediaObjects/"
                   "41586_2021_3493_MOESM6_ESM.xlsx")
LOCAL_XLSX = "data/validation/stukalov_supp_table2.xlsx"
SHEET = "A - Significant interactions"
PLK1_GENE = "PLK1"
TARGET_BAIT = "NSP13"


def download():
    urllib.request.urlretrieve(SUPP_TABLE_URL, LOCAL_XLSX)


def load_significant_interactions():
    wb = openpyxl.load_workbook(LOCAL_XLSX, read_only=True)
    ws = wb[SHEET]
    rows = ws.iter_rows(values_only=True)
    header = next(rows)
    return pd.DataFrame(rows, columns=header)


def main():
    print(f"Downloading Stukalov et al. 2021 Supplementary Table 2 from {SUPP_TABLE_URL} ...")
    download()
    df = load_significant_interactions()
    print(f"  total significant interaction rows (all viruses/baits): {len(df)}")

    sars2 = df[df["bait_organism"] == "SARS-CoV-2"]
    baits = sorted(sars2["bait_name"].unique())
    print(f"  SARS-CoV-2 baits profiled: {baits}")

    nsp13 = sars2[sars2["bait_name"] == TARGET_BAIT]
    nsp13_genes = sorted(nsp13["gene_name"].unique())
    print(f"\nSARS-CoV-2 {TARGET_BAIT} -- {len(nsp13_genes)} significant interactor(s) per Stukalov et al.: {nsp13_genes}")
    plk1_found = PLK1_GENE in nsp13_genes
    print(f"Is {PLK1_GENE} among them? {plk1_found}")

    # cross-virus context: the SARS-CoV (2003) NSP13 ortholog, profiled in the same study
    sars1_nsp13 = df[(df["bait_organism"] == "SARS-CoV") & (df["bait_name"] == TARGET_BAIT)]
    sars1_genes = sorted(sars1_nsp13["gene_name"].unique())
    print(f"\nSARS-CoV (2003) {TARGET_BAIT} ortholog -- {len(sars1_genes)} interactor(s): {sars1_genes}")
    print(f"Is {PLK1_GENE} among them? {PLK1_GENE in sars1_genes}")

    # context: how sparse is Nsp13 relative to other SARS-CoV-2 baits in this same dataset?
    counts = sars2.groupby("bait_name")["gene_name"].nunique().sort_values(ascending=False)
    print("\nInteractor count per SARS-CoV-2 bait in this independent dataset (context):")
    print(counts.to_string())

    out = pd.DataFrame([{
        "independent_source": "Stukalov et al. 2021 (Nature 594:246-252), Supplementary Table 2",
        "pmid": "33845483",
        "sars_cov2_nsp13_interactors": ";".join(nsp13_genes),
        "n_sars_cov2_nsp13_interactors": len(nsp13_genes),
        "plk1_reported_as_sars_cov2_nsp13_interactor": plk1_found,
        "sars_cov_2003_nsp13_interactors": ";".join(sars1_genes),
        "plk1_reported_as_sars_cov_2003_nsp13_interactor": PLK1_GENE in sars1_genes,
    }])
    out.to_csv("data/validation/nsp13_independent_crosscheck.csv", index=False)
    print("\nsaved data/validation/nsp13_independent_crosscheck.csv")


if __name__ == "__main__":
    main()
