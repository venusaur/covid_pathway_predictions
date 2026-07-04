# SARS-CoV-2 Host Interactome: Pathway & Missing-Component Analysis

Analysis built on the real protein-protein interaction (PPI) data from Gordon et al. 2020 (full citation and other data sources below).

## Data sources

All five are live public APIs queried at run time - no credentials required, but this means the scripts need internet access and will vary slightly if the underlying databases update.

| Source | Used for | Used by |
|---|---|---|
| [NDEx](https://www.ndexbio.org/) - network `HEK293T_SARS-CoV-2`, UUID `43803262-6d69-11ea-bfdc-0ac135e8bacf` | The real 332 SARS-CoV-2-human interactions (the dataset behind Fig. 3 of Gordon et al. 2020) | `fetch_network.py` |
| [Enrichr](https://maayanlab.cloud/Enrichr/) (via `gseapy`) - GO Biological Process, KEGG, Reactome libraries | Pathway enrichment per bait and overall; library gene-set coverage (denominator audit) | `enrichment.py`, `validation/02_coverage_audit.py`, `validation/03_strip_circular.py` |
| [STRING](https://string-db.org/) REST API | Network-expansion (`add_nodes`) guilt-by-association scoring to predict missing components | `predict_missing.py`, `validation/04_threshold_sweep.py`, `validation/05_null_model.py`, `validation/06_holdout_precision_recall.py` |
| [UniProt](https://www.uniprot.org/) REST API (`rest.uniprot.org`) | Identifier integrity check - confirming each prey's UniProt accession resolves to its recorded HGNC gene symbol | `validation/01_identifier_integrity.py` |
| [IntAct](https://www.ebi.ac.uk/intact/) via PSICQUIC (EBI) | Independent external validation - checking whether other SARS-CoV-2 interactome papers (Stukalov, Gao, Chen, Brandherm, Thorne, etc.) confirm a predicted edge | `validation/07_external_validation.py` |

The underlying study this whole analysis is built on:

> Gordon, D.E. et al. **"A SARS-CoV-2 protein interaction map reveals targets for drug repurposing."** *Nature* 583, 459-468 (2020). https://doi.org/10.1038/s41586-020-2286-9

The paper's Krogan-lab AP-MS experiments identified 332 high-confidence interactions between 27 SARS-CoV-2 viral proteins and human host proteins. This project pulls that real dataset (not simulated) and layers three things on top of it:

1. **Pathway enrichment** - which biological pathways/complexes each viral protein's host targets fall into.
2. **Missing-component prediction** - guilt-by-association network expansion to hypothesize additional host proteins that plausibly interact with each viral protein but weren't detected in the original screen.
3. **Adversarial validation** - 8 tests (identifier integrity, coverage audit, circularity, negative controls, threshold robustness, null model, held-out precision/recall, external confirmation) checking whether the predictions in (2) are real or artifacts.

## Pipeline

Run in order from the project root (with the venv activated):

```bash
source .venv/bin/activate

python3 fetch_network.py       # -> data/interactome.csv
python3 enrichment.py          # -> data/pathways_overall.csv, data/pathways_per_bait.csv
python3 predict_missing.py     # -> data/predicted_missing_components.csv
python3 visualize.py           # -> output/overall_top_pathways.png, output/bait_pathway_heatmap.png
python3 visualize_prediction.py  # -> output/predicted_nsp13_plk1.png, output/predicted_nsp4_tomm70.png
```

| Script | Purpose |
|---|---|
| `fetch_network.py` | Downloads the real NDEx network (CX format), parses nodes/edges, keeps only viral(bait)-human(prey) interactions, writes a flat interaction table with MiST/BFDR/AvgSpec confidence scores. |
| `enrichment.py` | Runs Enrichr pathway enrichment (GO BP, KEGG, Reactome) on the full 332-protein set and separately on each viral protein's own prey set (baits with < 3 preys are skipped). |
| `predict_missing.py` | For each bait with >= 5 known preys, queries STRING's network-expansion API (`add_nodes`) to find proteins most connected to that prey set. Any returned protein not already in the 332-interaction dataset is logged as a candidate "missing component," ranked by number of edges to known preys and STRING confidence. Also flags candidates that are already a confirmed interactor of a *different* viral bait, which is stronger independent evidence. |
| `visualize.py` | Bar chart of top overall enriched pathways; heatmap of bait x pathway significance (recreates the logic of the paper's Fig. 2b). |
| `visualize_prediction.py` | Small network diagrams for two illustrative predictions (see below). |

## Project layout

```
fetch_network.py
enrichment.py
predict_missing.py
visualize.py
visualize_prediction.py
validation/                        # 8-test adversarial validation suite (see below)
  01_identifier_integrity.py ... 07_external_validation.py
  VALIDATION_REPORT.md             # full per-test verdicts and honest caveats
data/                              # all CSVs (raw + computed tables)
  interactome.csv                  # 332 real viral-host interactions
  pathways_overall.csv             # enrichment across all host proteins
  pathways_per_bait.csv            # enrichment per viral protein
  predicted_missing_components.csv # guilt-by-association candidates per bait
  validation/                      # validation test outputs (CSVs)
output/                            # all plots
  overall_top_pathways.png
  bait_pathway_heatmap.png
  predicted_nsp13_plk1.png
  predicted_nsp4_tomm70.png
  validation/                      # validation test outputs (PNGs)
```

## Key findings (updated after adversarial validation -- see below)

**Pathway enrichment** recovers the paper's own biology as a sanity check - Reactome's own "SARS-CoV Infections" and "Potential Therapeutics For SARS" terms come up significant, alongside centrosome/mitotic pathways (driven by Nsp13), DNA replication (Nsp1), and mitochondrial protein import (Nsp4, Orf9b, Orf9c). But see the coverage audit below before trusting the ranking too literally.

**Missing-component predictions**, after 8 rounds of adversarial testing (`validation/`, full report in `validation/VALIDATION_REPORT.md`):

- **Nsp13 -> PLK1**: the stronger of the two calls. Not a generic STRING hub (0/25 hits under a degree-preserving null), stays specific under negative controls, and is stable across most of a reasonable confidence/width parameter band. Still **not independently confirmed** - no external SARS-CoV-2 interactome dataset (checked via IntAct) reports PLK1 binding Nsp13; PLK1's only independent SARS-CoV-2 record is with Spike.
- **Nsp4 -> TOMM70**: originally framed as "cross-validated" because TOMM70 is a real, confirmed Orf9b interactor. That framing **overstated it**. It passes the null model (not a generic hub) but fails the threshold-robustness sweep (only appears at the loosest STRING confidence setting - the exact fishing signature the test is designed to catch), and IntAct shows TOMM70 robustly confirmed across 7 independent papers **always with ORF9b, never with Nsp4**. Downgrade this to a weak, unconfirmed hypothesis, not a validated finding.

Quantitatively, the guilt-by-association method recovers a held-out true interaction 11.1% of the time at 3.5% precision on average, versus 0%/0% for a degree-matched null (real signal, above chance) - but recall varies from 0% (M, Nsp7, Nsp12, Orf9b) to 46.7% (N), so it works well for some baits and not at all for others.

## Caveat

Missing-component predictions are guilt-by-association on generic human PPI/co-expression data (STRING) - they indicate "this protein fits the same complex," not "this protein binds the virus." Treat the ranked list as hypotheses to prioritize for follow-up AP-MS/co-IP, the same way the original paper prioritized drug candidates rather than asserting mechanism outright. Every one of tests 1-7 in the validation suite still lives inside the STRING/Enrichr/UniProt ecosystem - only the external IntAct check (test 8) can actually confirm a prediction rather than merely fail to reject it, and it did not confirm either headline call.

## Validation suite

Two things the original analysis didn't do until challenged: quantify accuracy against a null baseline, and check independent data. Both live in `validation/`:

```bash
source .venv/bin/activate
python3 validation/01_identifier_integrity.py     # every gene symbol/UniProt accession checks out
python3 validation/02_coverage_audit.py           # Reactome only annotates 68.7% of our 332 genes
python3 validation/03_strip_circular.py           # ranking survives removing SARS-named Reactome terms
python3 validation/04_threshold_sweep.py          # Nsp13-PLK1 robust; Nsp4-TOMM70 fragile
python3 validation/05_null_model.py               # neither prediction is a generic hub artifact
python3 validation/06_holdout_precision_recall.py # 11.1%/3.5% recall/precision vs 0%/0% null
python3 validation/07_external_validation.py      # neither call independently confirmed via IntAct
```

Full writeup, per-test verdicts, and honest caveats about the validation's own limits: [`validation/VALIDATION_REPORT.md`](validation/VALIDATION_REPORT.md).
