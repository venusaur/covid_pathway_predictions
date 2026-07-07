"""Next-step 4: broaden pathway coverage before trusting the enrichment ranking.

Test 3's coverage audit found GO/KEGG/Reactome annotate only 68.7%/46.4%/91.6% of the
332 prey genes, unevenly across baits, and left 16 genes with zero annotation in any of
the three. This script adds three differently-curated libraries -- WikiPathways 2023,
BioPlanet 2019, and MSigDB Hallmark 2020 -- and asks:
  1. How many of the 16 previously-unannotated genes pick up any annotation now?
  2. Does overall and per-bait coverage even out, or do the same baits stay dark?
If the gaps persist across six independent libraries, they reflect genuinely
under-characterized biology, not one database's blind spot.
"""
import pandas as pd
import gseapy as gp

ORIGINAL_LIBS = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
ADDED_LIBS = ["WikiPathway_2023_Human", "BioPlanet_2019", "MSigDB_Hallmark_2020"]
ALL_LIBS = ORIGINAL_LIBS + ADDED_LIBS

PREVIOUSLY_UNANNOTATED = ['AAR2', 'ARL6IP6', 'C1orf50', 'CCDC86', 'CEP112', 'CEP350',
                          'CLCC1', 'CLIP4', 'FKBP15', 'MAP7D1', 'MIPOL1', 'PRRC2B',
                          'QSOX2', 'TMEM39B', 'USP54', 'ZC3H18']


def library_universe(name):
    gs = gp.get_library(name=name, organism="Human")
    u = set()
    for genes in gs.values():
        u.update(genes)
    return u


def main():
    gordon = pd.read_csv("data/interactome.csv")
    all_genes = set(gordon["prey_gene"].unique())

    universes = {}
    rows = []
    for lib in ALL_LIBS:
        print(f"Fetching {lib}...")
        uni = library_universe(lib)
        universes[lib] = uni
        covered = all_genes & uni
        rows.append({
            "library": lib,
            "added_in_this_step": lib in ADDED_LIBS,
            "our_genes_covered": len(covered),
            "coverage_fraction": round(len(covered) / len(all_genes), 3),
        })
    cov = pd.DataFrame(rows)
    cov.to_csv("data/validation/coverage_broadened.csv", index=False)
    print("\n=== Per-library coverage of the 332 prey genes ===")
    print(cov.to_string(index=False))

    # union coverage: original 3 vs all 6
    orig_union = set().union(*(universes[l] for l in ORIGINAL_LIBS))
    all_union = set().union(*(universes[l] for l in ALL_LIBS))
    orig_cov = len(all_genes & orig_union)
    all_cov = len(all_genes & all_union)
    print(f"\nUnion coverage, original 3 libraries: {orig_cov}/{len(all_genes)} ({orig_cov/len(all_genes):.1%})")
    print(f"Union coverage, all 6 libraries:      {all_cov}/{len(all_genes)} ({all_cov/len(all_genes):.1%})")

    # do the 16 previously-unannotated genes get rescued?
    print(f"\n=== The 16 previously-unannotated genes, checked against the 3 new libraries ===")
    rescued, still_dark = [], []
    rescue_rows = []
    for g in PREVIOUSLY_UNANNOTATED:
        hits = [lib for lib in ADDED_LIBS if g in universes[lib]]
        (rescued if hits else still_dark).append(g)
        rescue_rows.append({"gene": g, "rescued": bool(hits), "found_in": ";".join(hits)})
        print(f"  {g:10s} {'RESCUED by ' + ', '.join(hits) if hits else 'still unannotated'}")
    pd.DataFrame(rescue_rows).to_csv("data/validation/unannotated_rescue.csv", index=False)

    print(f"\nRescued by the new libraries: {len(rescued)}/16 -> {rescued}")
    print(f"Still unannotated across all 6 libraries: {len(still_dark)}/16 -> {still_dark}")

    # per-bait annotation depth (mean libraries-per-gene), all 6 libraries
    per_gene_libcount = {g: sum(1 for l in ALL_LIBS if g in universes[l]) for g in all_genes}
    gordon = gordon.assign(libcount=gordon["prey_gene"].map(per_gene_libcount))
    per_bait = gordon.groupby("bait")["libcount"].mean().sort_values()
    print("\n=== Per-bait mean annotation depth across all 6 libraries (lowest first) ===")
    print(per_bait.round(2).to_string())


if __name__ == "__main__":
    main()
