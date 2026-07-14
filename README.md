# SARS-CoV-2 Host Interactome: Pathway & Missing-Component Analysis

A computational reanalysis of the Gordon et al. 2020 SARS-CoV-2–human protein interactome, with
an unusually heavy emphasis on **not fooling yourself**: every prediction is put through a
14-test adversarial validation suite before it's believed.

> **Headline:** guilt-by-association network expansion looks validated under internal
> cross-validation (11.1% recall) but performs at chance against an independent dataset
> (0.07% vs 0.02% null). The apparent skill is literature-co-mention bias. One prediction
> (**Nsp13 → PLK1 / the centrosome**) survived every *association-based* test and reproduced on an
> independent physical-interaction network — but no bias-free line supports a *direct* interaction:
> no docking motif, not a CRISPR host factor, no co-IP/imaging in four years of literature. The AF3
> structural run (11 jobs, completed 2026-07-12) is **decisive**: with its pipeline control validated
> (Nsp10:Nsp16, ipTM 0.79), it excludes a direct interface for **every** candidate — PLK1, CEP135,
> PRKACA, AKAP9, and the genome-wide SPRINT top hits MCM7/TRAF2 (all ≤ 0.25 vs 0.79). So **no direct
> human partner for Nsp13 is supported**; the co-IP-validated **PKA (PRKACA)** contact is best read as
> **indirect/bridged**, and the SPRINT ranks are confirmed as out-of-distribution homology artifacts.
> The one open lead is a *kinase–substrate* relationship (4 PLK1-consensus sites in Nsp13), a bench
> assay. A second prediction (Nsp4 → TOMM70) was **withdrawn** on scrutiny.

---

## Reading guide

The narrative lives in a few documents; read them in this order depending on how deep you want to go.

| Read | What it is | For |
|---|---|---|
| **[`report/HOW_WE_GOT_HERE.md`](report/HOW_WE_GOT_HERE.md)** | **Plain-language walkthrough** — every step in the order we did it, the method behind each, and why it mattered. No virology background needed. | Anyone (easiest start) |
| **[`report/NSP13_CENTROSOME_SUMMARY.md`](report/NSP13_CENTROSOME_SUMMARY.md)** | **Single-file synthesis of the whole investigation** — network → SPRINT → AF3 (excludes direct binding for all candidates) → PKA relocated to indirect → wet-lab. The fastest complete read. | The technical synthesis |
| **[`report/MANUSCRIPT_DRAFT.md`](report/MANUSCRIPT_DRAFT.md)** | The whole study as a **publishable negative-result/methods paper**: abstract, methods, results, discussion, limitations, references. | Writing it up |
| **[`report/q1_prediction_validation_rate.md`](report/q1_prediction_validation_rate.md)** | Interactome-scale companion: what fraction of the 279 PIPE4/SPRINT SARS-CoV-2 predictions validate experimentally (28.3%, 1.53× over null), and the patterns in what holds up | Broader context on prediction yield |
| **[`report/VALIDATION_REPORT.md`](report/VALIDATION_REPORT.md)** | The 14-test adversarial suite — per-test verdicts and caveats on the validation's own limits | The network methodology in detail |
| **[`report/FUNCTIONAL_GENOMICS.md`](report/FUNCTIONAL_GENOMICS.md)** · **[`report/EXPERIMENTAL_DESIGN.md`](report/EXPERIMENTAL_DESIGN.md)** | The CRISPR host-factor line, and the PKA/AKAP9-led wet-lab co-IP design | Phenotype line + what's next |
| **This README** | Overview, repo map, how to run, condensed findings | Orientation + reproduction |

All analysis reports and writeups live in [`report/`](report/). The `experimental/` folder holds the remaining working files:
- [`alphafold/ALPHAFOLD_PROTOCOL.md`](experimental/alphafold/ALPHAFOLD_PROTOCOL.md): the structural run matrix, interpretation framework, and **the 2026-07-12 results — with the pipeline control validated (0.79), AF3 excludes a direct interface for all 8 candidates (PLK1, CEP135, PRKACA, AKAP9, MCM7, TRAF2); PKA relocated to indirect**
- [`alphafold/ALPHAFOLD_TUTORIAL.md`](experimental/alphafold/ALPHAFOLD_TUTORIAL.md): hands-on, run-it-yourself AF3 walkthrough (no experience needed)
- [`alphafold/alphafold_inputs/scorecard.py`](experimental/alphafold/alphafold_inputs/scorecard.py): scores the AF3 output zips against the controls and prints the ipTM scorecard
- [`plk1_motif_scan.py`](experimental/plk1_motif_scan.py): scans Nsp13 for the PLK1 Polo-box docking + kinase-substrate motifs
- [`literature/`](experimental/literature/): published-literature dive against the open questions, with a full source log ([`LITERATURE_REVIEW.md`](experimental/literature/LITERATURE_REVIEW.md), [`SOURCE_LIST.md`](experimental/literature/SOURCE_LIST.md), [`SOURCES.md`](experimental/literature/SOURCES.md))

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
validation/             the 14-test adversarial suite scripts (01..13 + inline negative controls)
experimental/           post-ceiling follow-up: AlphaFold, motif scan, literature
  alphafold/              all AlphaFold3 work:
    ALPHAFOLD_PROTOCOL.md   run matrix, interpretation, results
    ALPHAFOLD_TUTORIAL.md   hands-on run-it-yourself guide
    alphafold_inputs/       FASTAs + af3_server_entries.txt + scorecard.py
    alphafold_outputs/      AF3 result zips + scored metrics
  plk1_motif_scan.py      Nsp13 PLK1-motif scan
  literature/             published-literature dive + source logs
data/                   inputs + all computed tables (incl. q1_validation/)
  validation/             per-test CSV outputs
output/                 all figures (PNG)
  validation/             per-test figures
report/                 all analysis reports & writeups:
  HOW_WE_GOT_HERE.md, NSP13_CENTROSOME_SUMMARY.md, MANUSCRIPT_DRAFT.md,
  VALIDATION_REPORT.md, sprint_nsp13_plk1_validation.md,
  q1_prediction_validation_rate.md, FUNCTIONAL_GENOMICS.md,
  EXPERIMENTAL_DESIGN.md, SPRINT_DATA.md
SPRINT_DATA/            copied SPRINT scan results + README
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
- **Nsp13 → PLK1 — direct interaction unsupported; redirected to PKA/AKAP9.** It survived every
  *association-based* test and reproduced on IntAct's independent physical-interaction graph
  (6/40 Nsp13 preys are direct PLK1 partners, p=0.0005) — not a STRING/text-mining artifact. But no
  bias-free line supports a *direct* contact: no PLK1 docking motif in Nsp13, not a CRISPR
  host-factor hit, and no co-IP/imaging in four years places Nsp13 at PLK1 or the centrosome. The
  AF3 structural run meant to settle it is **decisive** — with its pipeline control validated
  (Nsp10:Nsp16, 0.79), it scores every Nsp13 pair (PLK1×3, CEP135, PRKACA, AKAP9, and the genome-wide
  SPRINT top hits MCM7/TRAF2) at or below the 0.21 decoy floor. The old "AF3 is blind to Nsp13" read
  does not survive that control: CEP135/PRKACA are co-IP/association evidence, which does not prove a
  *direct* interface, so their nulls are expected for indirect partners. Net: **no direct human
  partner for Nsp13 is supported.** The one route to the centrosome with an *experimentally validated*
  Nsp13 contact is **PKA** (Nsp13↔PRKACA, Wang et al., PMC11019953; PKA anchored via AKAP9, a Gordon
  prey) — but AF3's null relocates it to an **indirect/bridged** association, directness now an open
  question. **The wet-lab experiments are now decisive**: a direct-vs-bridged co-IP/reconstitution for
  PKA, and an in-vitro kinase assay on Nsp13's four PLK1-consensus sites (the one open lead).
- **Nsp4 → TOMM70**: **withdrawn.** Threshold-fragile, TOMM70's real viral partner is ORF9b (7
  papers), and it fails to reproduce on IntAct. Four independent lines agree.

**Caveat.** These are guilt-by-association predictions: "fits the same complex," not "binds the
virus." Treat them as hypotheses to prioritize, exactly as the original paper treated its drug
candidates. Only checks against independent data can *confirm* rather than fail-to-reject;
none did.

**Do the predictions hold up at scale?** Beyond the two headline cases, we audited a whole published
PIPE4/SPRINT SARS-CoV-2 predicted interactome (279 high-confidence pairs) against experimental
evidence: **28.3% are experimentally supported, 1.53× over a viral-protein-matched null (p<0.0005)**,
real but modest. SPRINT's confidence tracks validation (19% to 39% across score tertiles) while
PIPE4's is flat, so SPRINT-ranked or two-method-consensus filtering is the practical lever. Full
audit: [`report/q1_prediction_validation_rate.md`](report/q1_prediction_validation_rate.md).

---

## Future direction

The computational program is complete: no direct human partner for Nsp13 is supported, and the
interactome-scale audit shows why convergent predictors mislead on viral baits (shared homology and
hub signal, not interfaces). Three tracks remain, in priority order:

1. **Publish the negative.** `report/MANUSCRIPT_DRAFT.md` is a near-complete negative-result/methods
   paper; the remaining gaps are an abstract and turning the figure legends into panels.
2. **Take the one live lead to the bench.** Nsp13 carries four PLK1 kinase-substrate consensus sites
   (T144, T199, S377, S385), a *kinase–substrate* relationship AF3 cannot see; an in-vitro PLK1
   kinase assay tests it. Separately, a direct-vs-bridged co-IP and reconstitution settles whether the
   co-IP-validated PKA contact is direct (`report/EXPERIMENTAL_DESIGN.md`).
3. **Reuse the machinery on a new target.** The validated-control AF3 scorecard, the genome-wide
   SPRINT scan, and the validation-rate audit transfer to any other viral protein with no rework.

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
