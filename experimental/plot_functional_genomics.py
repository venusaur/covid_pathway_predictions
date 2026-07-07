"""Figure: PLK1 vs ACE2 across three genome-wide SARS-CoV-2 CRISPR KO screens.
Reads the cached supplementary tables (fetched by functional_genomics_check.py)."""
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def rank_table(fp, sheet, rank_col):
    ws = openpyxl.load_workbook(fp, read_only=True)[sheet]
    rows = [r for r in list(ws.iter_rows(values_only=True))[1:] if r[0]]
    return {r[0]: r[rank_col] for r in rows}, len(rows)


def pct(rank, n):
    return 100 * (1 - (rank - 1) / n)


def main():
    calu, ncalu = rank_table("/tmp/fg_calu3.xlsx", "Calu3_Gattinara_avg_zscore", 2)
    vero, nvero = rank_table("/tmp/fg_vero.xlsx", "VeroE6_avg_zscore", 2)
    bier, nb = rank_table("/tmp/fg_biering.xlsx", "Supplementary Table 1", 3)

    screens = [("Calu-3 (lung)", calu, ncalu), ("Vero E6", vero, nvero),
               ("Biering bidirectional", bier, nb)]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    for y, (name, d, n) in enumerate(screens):
        ax.barh(y, 1.0, left=99, height=0.5, color="#e7c9c4", zorder=1)  # top-1% hit zone
        ax.plot([0, 100], [y, y], color="#ccc", lw=1, zorder=0)
        p_ace, p_plk = pct(d["ACE2"], n), pct(d["PLK1"], n)
        ax.scatter(p_ace, y, s=120, color="#3b6fa0", zorder=3,
                   label="ACE2 (top host factor)" if y == 0 else "")
        ax.scatter(p_plk, y, s=140, color="#c0392b", marker="D", zorder=3,
                   label="PLK1" if y == 0 else "")
        ax.annotate(f"#{d['ACE2']}", (p_ace, y), xytext=(0, 9), textcoords="offset points",
                    ha="center", fontsize=8, color="#3b6fa0")
        ax.annotate(f"#{d['PLK1']}", (p_plk, y), xytext=(0, -16), textcoords="offset points",
                    ha="center", fontsize=8, color="#c0392b")

    ax.set_yticks(range(len(screens)))
    ax.set_yticklabels([s[0] for s in screens])
    ax.set_xlim(0, 102)
    ax.set_xlabel("Pro-viral enrichment percentile  (100 = strongest host-dependency factor)")
    ax.set_title("PLK1 is not a host-factor hit in any genome-wide CRISPR screen\n"
                 "ACE2 sits at the ceiling in all three; PLK1 never enters the top-1% hit zone (shaded)")
    ax.text(99.5, len(screens) - 0.4, "top-1%\nhit zone", ha="center", va="top",
            fontsize=7, color="#b23a2e")
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    ax.invert_yaxis()
    ax.grid(True, axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig("output/validation/functional_genomics_plk1.png", dpi=150)
    print("saved output/validation/functional_genomics_plk1.png")
    for name, d, n in screens:
        print(f"  {name}: PLK1 percentile {pct(d['PLK1'], n):.1f}")


if __name__ == "__main__":
    main()
