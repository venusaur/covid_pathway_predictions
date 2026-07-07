"""Pre-flight data-integrity gate for the SARS-CoV-2 interactome pipeline.

Run this BEFORE enrichment.py / predict_missing.py. It fails loudly (exit code 1) on the
two classes of silent corruption the validation suite had to catch after the fact:
  - broken joins / malformed identifiers in the interactome (validation test 1), and
  - near-total loss of pathway-library coverage (validation test 3),
so any future analysis built on this pipeline gets those checks by default.

Usage:
    python3 preflight.py            # structural + coverage gates
    python3 preflight.py --offline  # structural gate only (no network)

    python3 preflight.py && python3 enrichment.py   # gate the pipeline
"""
import sys
import re
import pandas as pd

INTERACTOME = "data/interactome.csv"
REQUIRED_COLS = ["bait", "prey_gene", "prey_uniprot", "MIST", "BFDR", "AvgSpec"]
UNIPROT_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
COVERAGE_LIBS = ["GO_Biological_Process_2023", "Reactome_2022"]
COVERAGE_FLOOR = 0.30  # fail if a library annotates < 30% of prey genes (near-total loss)

GREEN, RED, YELLOW, RESET = "\033[32m", "\033[31m", "\033[33m", "\033[0m"


def fail(msg):
    print(f"{RED}PRE-FLIGHT FAIL:{RESET} {msg}")
    sys.exit(1)


def ok(msg):
    print(f"{GREEN}  ok{RESET}  {msg}")


def warn(msg):
    print(f"{YELLOW}  warn{RESET} {msg}")


def structural_gate():
    print("Structural gate (offline):")
    try:
        df = pd.read_csv(INTERACTOME)
    except FileNotFoundError:
        fail(f"{INTERACTOME} not found -- run fetch_network.py first.")

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        fail(f"interactome missing required columns: {missing}")
    ok(f"all required columns present ({len(df)} rows)")

    if df["prey_gene"].isna().any():
        fail(f"{int(df['prey_gene'].isna().sum())} rows have a null prey_gene (broken join).")
    if df["prey_uniprot"].isna().any():
        fail(f"{int(df['prey_uniprot'].isna().sum())} rows have a null prey_uniprot (broken join).")
    ok("no null gene / accession values")

    dups = int(df.duplicated(subset=["bait", "prey_gene"]).sum())
    if dups:
        fail(f"{dups} duplicate (bait, prey_gene) rows -- likely a bad merge.")
    ok("no duplicate (bait, prey) rows")

    bad = df.loc[~df["prey_uniprot"].str.match(UNIPROT_RE), "prey_uniprot"].unique()
    if len(bad):
        fail(f"{len(bad)} malformed UniProt accessions, e.g. {list(bad)[:5]} -- "
             f"a mangled join won't announce itself.")
    ok("all UniProt accessions well-formed")
    return df


def coverage_gate(df):
    print("Coverage gate (online):")
    try:
        import gseapy as gp
    except ImportError:
        warn("gseapy not installed -- skipping coverage gate.")
        return
    genes = set(df["prey_gene"].unique())
    for lib in COVERAGE_LIBS:
        try:
            gs = gp.get_library(name=lib, organism="Human")
        except Exception as ex:
            warn(f"could not fetch {lib} ({ex}) -- skipping.")
            continue
        universe = set().union(*gs.values())
        frac = len(genes & universe) / len(genes)
        if frac < COVERAGE_FLOOR:
            fail(f"{lib} annotates only {frac:.1%} of prey genes (< {COVERAGE_FLOOR:.0%} floor) "
                 f"-- enrichment would describe the library, not the biology.")
        ok(f"{lib}: {frac:.1%} coverage")


def main():
    offline = "--offline" in sys.argv
    print("=== SARS-CoV-2 interactome pipeline pre-flight ===")
    df = structural_gate()
    if not offline:
        coverage_gate(df)
    print(f"{GREEN}PRE-FLIGHT PASS{RESET} -- safe to run enrichment / prediction.")


if __name__ == "__main__":
    main()
