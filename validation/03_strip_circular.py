"""Test 5: remove Reactome/GO terms that were curated FROM the SARS-CoV-2 literature
itself (i.e. terms whose name references SARS/COVID/Coronavirus) and re-rank.
If the top hits survive without them, the result stands on independent signal.
If the ranking scrambles, some of the 'validation' was the paper being quoted
back through a curation layer.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CIRCULAR_PATTERN = r"SARS|COVID|Coronavirus"


def strip(df):
    mask = df["Term"].str.contains(CIRCULAR_PATTERN, case=False, regex=True)
    print(f"Removing {mask.sum()} circular terms out of {len(df)}:")
    print(df.loc[mask, "Term"].to_string(index=False))
    return df[~mask].copy()


def main():
    overall = pd.read_csv("data/pathways_overall.csv")
    print("=== BEFORE stripping: top 10 overall ===")
    print(overall.head(10)[["Term", "Adjusted P-value"]].to_string(index=False))

    stripped = strip(overall)
    stripped.to_csv("data/validation/pathways_overall_stripped.csv", index=False)

    print("\n=== AFTER stripping: top 10 overall ===")
    print(stripped.head(10)[["Term", "Adjusted P-value"]].to_string(index=False))

    # did the previously-#1 non-circular term move?
    before_ranks = {t: i for i, t in enumerate(overall["Term"])}
    after_ranks = {t: i for i, t in enumerate(stripped["Term"])}
    shifts = []
    for t in stripped["Term"].head(15):
        shifts.append((t, before_ranks[t], after_ranks[t]))
    print("\n=== Rank shift for terms now in the new top 15 (before_rank -> after_rank) ===")
    for t, b, a in shifts:
        print(f"  {b:>4} -> {a:>4}   {t}")

    fig, ax = plt.subplots(figsize=(9, 6))
    top = stripped.head(15).iloc[::-1]
    vals = -np.log10(top["Adjusted P-value"])
    ax.barh(top["Term"].str.slice(0, 55), vals, color="#3b6fa0")
    ax.set_xlabel("-log10(adjusted P-value)")
    ax.set_title("Top pathways with SARS/COVID-curated Reactome terms REMOVED\n"
                  "(tests whether ranking survives without circular annotation)")
    fig.tight_layout()
    fig.savefig("output/validation/overall_top_pathways_stripped.png", dpi=150)
    print("\nsaved output/validation/overall_top_pathways_stripped.png")


if __name__ == "__main__":
    main()
