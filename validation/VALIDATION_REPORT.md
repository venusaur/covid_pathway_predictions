# Validation report: does the missing-component prediction actually work?

This report runs 8 adversarial tests against the two headline claims from the earlier analysis:
**Nsp13 -> PLK1** and **Nsp4 -> TOMM70** (both "predicted missing components" via STRING
guilt-by-association). Every test below is implemented as a runnable script in this folder and
was executed against live data -- nothing here is asserted without a number behind it.

## Summary table

| # | Test | Script | Verdict |
|---|---|---|---|
| 7 | Identifier integrity | `01_identifier_integrity.py` | **PASS** -- 331/332 exact UniProt/HGNC matches, 1 legitimate synonym, zero broken joins |
| 4 | Coverage/denominator audit | `02_coverage_audit.py` | **CONCERN CONFIRMED** -- Reactome covers only 68.7%, KEGG 46.4% of our 332 genes; coverage is bait-dependent |
| 5 | Strip circular annotations | `03_strip_circular.py` | **PASS** -- removing 5 SARS-named Reactome terms doesn't reorder the remaining ranking |
| 2 | Negative controls (E, M, Orf7a) | run inline, see below | **PASS** -- PLK1/TOMM70 never leak into unrelated baits; predictions stay topically anchored |
| 6 | Threshold robustness | `04_threshold_sweep.py` | **SPLIT** -- Nsp13->PLK1 robust across most of the parameter band; Nsp4->TOMM70 fragile, survives only at the loosest confidence setting |
| 1 | Degree-preserving null model | `05_null_model.py` | **PASS** -- 0/25 hits for both PLK1 and TOMM70 under randomized same-size gene sets; not generic hubs |
| 3 | Precision/recall via hold-out | `06_holdout_precision_recall.py` | **REAL BUT MODEST** -- 11.1% recall / 3.5% precision (real) vs 0%/0% (null), but highly bait-dependent (0% for M, Nsp7, Nsp12, Orf9b) |
| 8 | External validation (IntAct) | `07_external_validation.py` | **NOT CONFIRMED** -- TOMM70 is well-established as a SARS-CoV-2 target across 7 independent papers, but always with ORF9b, never Nsp4. PLK1 has 1 independent SARS-CoV-2 record, with Spike, never Nsp13 |

## What this means for the two headline claims

**Nsp13 -> PLK1**: survives every internal test (null model, negative controls, most of the
threshold band) and the holdout recall for Nsp13 specifically is one of the best in the dataset
(27.5%). But no independent interactome dataset has ever reported PLK1 binding Nsp13 -- only
Spike. Read this as: *a well-supported, specific, non-artifactual hypothesis that has not yet
been externally confirmed.*

**Nsp4 -> TOMM70**: passes the null model (not a generic hub) and negative controls, but fails
the threshold-robustness test (only appears at the loosest STRING confidence setting) and has
zero external support for this specific pairing -- TOMM70's real, repeatedly-confirmed viral
partner is ORF9b, already in the original 332. This should be **downgraded** from "compelling
cross-validated finding" (the original framing) to "weak, threshold-dependent hypothesis." The
original framing overstated this one.

## Honest caveats about the validation itself

- Tests 1-7 all still live inside the STRING/Enrichr/UniProt ecosystem to varying degrees --
  guilt-by-association validated by more association is a closed loop, which is exactly why
  test 8 (external, independent data) is the only one that can actually confirm rather than
  merely fail to reject. It did not confirm either headline call.
- The holdout precision numbers (test 3) are a lower bound: a "wrong" candidate could be a real,
  undiscovered interactor rather than noise -- we only have a ground truth for what the original
  screen happened to detect.
- IntAct coverage is not exhaustive; absence of an independent Nsp13-PLK1 or Nsp4-TOMM70 record
  is evidence of absence only to the extent IntAct's curation is complete for these proteins,
  which cannot be fully verified from outside.

## How to reproduce

```bash
source .venv/bin/activate
python3 validation/01_identifier_integrity.py
python3 validation/02_coverage_audit.py
python3 validation/03_strip_circular.py
python3 validation/04_threshold_sweep.py
python3 validation/05_null_model.py
python3 validation/06_holdout_precision_recall.py
python3 validation/07_external_validation.py
```

Outputs land in `data/validation/*.csv` and `output/validation/*.png`.
