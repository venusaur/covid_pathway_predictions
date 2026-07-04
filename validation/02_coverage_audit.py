"""Test 4 (coverage / denominator audit): before trusting any enrichment bar chart,
find out how many of our 332 prey genes are even annotated in the reference
libraries we queried. If GO/Reactome only annotate a fraction of them, the
enrichment landscape partly reflects the library's coverage, not the biology.
"""
import pandas as pd
import gseapy as gp

LIBRARIES = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]


def main():
    df = pd.read_csv("data/interactome.csv")
    all_genes = set(df["prey_gene"].unique())
    print(f"Total unique prey genes: {len(all_genes)}")

    rows = []
    per_gene_hits = {g: set() for g in all_genes}
    for lib in LIBRARIES:
        print(f"Fetching library '{lib}'...")
        gene_sets = gp.get_library(name=lib, organism="Human")
        annotated_universe = set()
        for term, genes in gene_sets.items():
            annotated_universe.update(genes)
        covered = all_genes & annotated_universe
        for g in covered:
            per_gene_hits[g].add(lib)
        rows.append({
            "library": lib,
            "n_terms": len(gene_sets),
            "library_universe_size": len(annotated_universe),
            "our_genes_covered": len(covered),
            "our_genes_total": len(all_genes),
            "coverage_fraction": round(len(covered) / len(all_genes), 3),
        })

    cov_df = pd.DataFrame(rows)
    cov_df.to_csv("data/validation/coverage_audit.csv", index=False)
    print("\n=== Per-library coverage of our 332 prey genes ===")
    print(cov_df.to_string(index=False))

    any_lib = sum(1 for g, libs in per_gene_hits.items() if libs)
    all_libs = sum(1 for g, libs in per_gene_hits.items() if len(libs) == len(LIBRARIES))
    none = sorted(g for g, libs in per_gene_hits.items() if not libs)
    print(f"\nGenes annotated in at least one library: {any_lib}/{len(all_genes)} "
          f"({any_lib/len(all_genes):.1%})")
    print(f"Genes annotated in ALL three libraries: {all_libs}/{len(all_genes)} "
          f"({all_libs/len(all_genes):.1%})")
    print(f"Genes annotated in NONE of the three libraries ({len(none)}): {none}")

    per_gene_df = pd.DataFrame([
        {"gene": g, "n_libraries_annotated": len(libs), "libraries": ";".join(sorted(libs))}
        for g, libs in per_gene_hits.items()
    ]).sort_values("n_libraries_annotated")
    per_gene_df.to_csv("data/validation/coverage_per_gene.csv", index=False)


if __name__ == "__main__":
    main()
