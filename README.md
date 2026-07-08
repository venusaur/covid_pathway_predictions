# SARS-CoV-2 Host Interactome: Pathway & Missing-Component Analysis

A computational reanalysis of the Gordon et al. 2020 SARS-CoV-2–human protein interactome, with
an unusually heavy emphasis on **not fooling yourself**: every prediction is put through a
14-test adversarial validation suite before it's believed.

> **Headline:** guilt-by-association network expansion looks validated under internal
> cross-validation (11.1% recall) but performs at chance against an independent dataset
> (0.07% vs 0.02% null). The apparent skill is literature-co-mention bias. One prediction
> survived everything computation can throw at it (**Nsp13 → PLK1 / the centrosome**), reproducing
> even on an independent physical-interaction network, but it hits a hard confirmation ceiling
> and now needs a wet-lab co-IP. A second (Nsp4 → TOMM70) was **withdrawn** on scrutiny.

---

## Repository structure

```
README.md               ← you are here
requirements.txt
pipeline/               core analysis, run in order from the project root
  fetch_network.py        NDEx  -> data/interactome.csv (the real 332 interactions)
  enrichment.py           Enrichr pathway enrichment -> data/pathways_*.csv
  predict_missing.py      STRING guilt-by-association -> data/predicted_missing_components.csv
  visualize.py            overall + per-bait pathway figures
  visualize_prediction.py network diagrams for the two headline predictions
  preflight.py            data-integrity gate (run before enrichment/prediction)
validation/             the 14-test adversarial suite (01..13 + inline negative controls)
  VALIDATION_REPORT.md    full per-test writeup
experimental/           post-ceiling follow-up: co-IP design, AlphaFold, functional genomics
  alphafold_inputs/       ready-to-run FASTAs + query_sequences.txt
data/                   inputs + all computed tables
  validation/             per-test CSV outputs
output/                 all figures (PNG)
  validation/             per-test figures
report/                 the polished findings report (HTML + PDF)
```

Everything runs against **live public APIs** (no credentials), so scripts need internet access
and results can shift slightly if the underlying databases update. Two large re-downloadable
intermediates (`data/huri_edges.csv`, the Stukalov supplement) are git-ignored and regenerated on
demand.

---

## Quickstart

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# core pipeline, run from the project root, in order
python3 pipeline/fetch_network.py       # -> data/interactome.csv
python3 pipeline/preflight.py           # data-integrity gate (fails loudly on a broken join)
python3 pipeline/enrichment.py          # -> data/pathways_overall.csv, pathways_per_bait.csv
python3 pipeline/predict_missing.py     # -> data/predicted_missing_components.csv
python3 pipeline/visualize.py           # -> output/*.png
python3 pipeline/visualize_prediction.py

# validation suite (each writes to data/validation/ and output/validation/)
for f in validation/[0-9]*.py; do python3 "$f"; done
```

All scripts assume the **project root** as the working directory.

---

## Findings (condensed; full version in the report)

**Pathway enrichment** recovers the paper's own biology (centrosome/mitosis on Nsp13, DNA
replication on Nsp1, mitochondrial import on Nsp4/Orf9b/Orf9c). But GO/KEGG/Reactome annotate only
68.7%/46.4%/91.6% of the 332 genes, unevenly across baits, and broadening to six libraries barely
helps (95.2%→96.4%; 12 of 16 unannotated genes stay dark). Cross-bait comparisons are never quite
on equal footing.

**The method mostly doesn't generalize.** Internal holdout suggested skill (11.1% recall vs
0% null), but against Stukalov et al.'s independent AP-MS as ground truth, mean recall collapsed to
0.07% (vs 0.02% null) across 18 baits: chance. The diagnosis: holdout recall tracks how densely a
bait's preys are already clustered in STRING (ρ +0.65), and that clustering is driven by the
literature *text-mining* channel (ρ +0.78) more than hard experimental evidence (ρ +0.50). The
method rewards fame, not evidence, so it looks good on already-famous proteins and fails on novel
biology.

**The two headline predictions:**
- **Nsp13 → PLK1**: survived every internal test *and* reproduced on IntAct's independent
  physical-interaction graph (6/40 Nsp13 preys are direct PLK1 partners, p=0.0005). Not a
  STRING/text-mining artifact. But it's not directly confirmed in any viral interactome, HuRI
  (unbiased Y2H) is blind to the centrosome so can't adjudicate it, and it's not a CRISPR
  host-factor hit, so the biology, if real, is centrosome/cell-cycle perturbation, not an
  essential host factor. Motif analysis (no PLK1 docking motif) points to a **bridged model**
  (Nsp13 → CEP scaffolds → PLK1). Confirmation now requires a wet-lab co-IP.
- **Nsp4 → TOMM70**: **withdrawn.** Threshold-fragile, TOMM70's real viral partner is ORF9b (7
  papers), and it fails to reproduce on IntAct. Four independent lines agree.

**Caveat.** These are guilt-by-association predictions: "fits the same complex," not "binds the
virus." Treat them as hypotheses to prioritize, exactly as the original paper treated its drug
candidates. Only checks against independent data can *confirm* rather than fail-to-reject;
none did.

---

## Data sources

All are live public APIs / repositories queried at run time.

| Source | Used for |
|---|---|
| [NDEx](https://www.ndexbio.org/), `HEK293T_SARS-CoV-2` (`43803262-…`) and HuRI (`73bc2c06-…`) | The real 332 interactions; the HuRI Y2H graph |
| [Enrichr](https://maayanlab.cloud/Enrichr/) (via `gseapy`) | Pathway enrichment + coverage audits |
| [STRING](https://string-db.org/) REST API | Guilt-by-association scoring, null model, thresholds, evidence-channel diagnosis |
| [UniProt](https://www.uniprot.org/) REST API | Identifier integrity; sequences for the experimental work |
| [IntAct](https://www.ebi.ac.uk/intact/) via PSICQUIC | Independent external validation; the alternative physical-interaction substrate |
| [Springer/Nature](https://www.nature.com/) static content | Independent AP-MS (Stukalov 2021) and CRISPR screen supplements |

Underlying study:

> Gordon, D.E. et al. **"A SARS-CoV-2 protein interaction map reveals targets for drug repurposing."**
> *Nature* 583, 459–468 (2020). https://doi.org/10.1038/s41586-020-2286-9
