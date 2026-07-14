# covid_ppi_pathways

A computational reanalysis of the SARS-CoV-2 to human protein interactome. It started as a
guilt-by-association pipeline to expand the Gordon et al. 2020 AP-MS map, and turned into an
adversarial audit of what those predictions are worth.

## Bottom line

**Nsp13 has no direct, high-confidence human protein partner among any candidate we raised.**

One prediction, **Nsp13 to PLK1 / the centrosome**, survived every association-based test
(STRING+IntAct network, p=0.0005; SPRINT sequence scan). But every bias-free line declined to
support a *direct* interaction:

- **Structure (AlphaFold3):** with a validated pipeline control (Nsp10:Nsp16 folded at ipTM 0.79),
  AF3 placed PLK1, CEP135, PRKACA, AKAP9, and the genome-wide SPRINT top hits (MCM7, TRAF2) at or
  below the decoy floor. A direct interface is excluded for all of them, not just PLK1.
- **Phenotype (CRISPR):** PLK1 is not a host-dependency factor in three genome-wide knockout screens.
- **Literature (2020 to 2025):** no co-IP or imaging study places Nsp13 at PLK1 or the centrosome.

The co-IP-validated **PKA (PRKACA)** contact is best modelled as **indirect / bridged**, consistent
with its own AF3 null. The one open computational lead is not a binding partner but a
**kinase-substrate** relationship (4 PLK1-consensus Ser/Thr sites in Nsp13), which is a bench assay.

Alongside the negative, we audited the source predictions at scale: of 279 high-confidence
PIPE4/SPRINT SARS-CoV-2 to human predictions, **28.3% (79/279) are backed by an experimental
BioGRID interaction**, a 1.53x enrichment over a viral-matched permutation null (18.6%, p<0.0005).
SPRINT scores are calibrated (validation rises with score tertile: 19/27/39%); PIPE4 scores are
flat. Validation concentrates on abundant compartments (mitochondrion 64%, ribosome 47%), which is
partly an AP-MS/BioID background signal.

## Repo layout

```
pipeline/       Prediction pipeline: fetch network, enrichment, predict missing components, visualize
validation/     14-test adversarial validation suite (01..13, numbered by run order)
experimental/   Follow-up lines that escape the network's biases:
  alphafold/      AF3 docking protocol, inputs, outputs, scorecard
  literature/     Manual literature review
  *_motif_scan, functional_genomics_check   sequence + CRISPR-phenotype checks
data/           Interactome tables, HuRI edges, and data/q1_validation/ (the 279-prediction audit)
SPRINT_DATA/    SPRINT sequence-homology scan inputs, scans, results
output/         Figures (heatmaps, pathway plots, predicted-pair visualizations)
report/         All write-ups (start here)
```

## Reports (start here)

| Document | What it is |
|---|---|
| `report/NSP13_CENTROSOME_SUMMARY.md` | Single-file synthesis of the whole investigation. Read first. |
| `report/HOW_WE_GOT_HERE.md` | Plain-language, step-by-step walkthrough of how we got the result. |
| `report/MANUSCRIPT_DRAFT.md` | Manuscript draft with figures and methods. |
| `report/q1_prediction_validation_rate.md` | The 279-prediction validation-rate audit. |
| `report/VALIDATION_REPORT.md` | The 14-test adversarial validation of the network predictions. |
| `report/sprint_nsp13_plk1_validation.md` | SPRINT sequence scan + PLK1 motif addendum. |
| `report/FUNCTIONAL_GENOMICS.md` | CRISPR host-factor triangulation. |
| `report/EXPERIMENTAL_DESIGN.md` | Wet-lab plan (co-IP, direct-vs-bridged, imaging). |

## Setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Most scripts are stdlib-only (numpy/scipy where noted). To reproduce the Q1 audit:

```bash
python data/q1_validation/q1_analysis.py          # reproduces all Q1 numbers
python data/q1_validation/q1_analysis.py --test    # self-check
```

## Future direction

Three tracks, in priority order:

1. **Publish the negative.** `MANUSCRIPT_DRAFT.md` is close; it still needs an abstract and the
   primary Borealis citation. The methodological point (network guilt-by-association looks skilled
   under internal cross-validation but collapses against an independent dataset) stands on its own.
2. **Take the surviving leads to the bench.** The 4 PLK1-consensus kinase-substrate sites in Nsp13,
   and the direct-vs-bridged question for the PKA/AKAP9 contact. `EXPERIMENTAL_DESIGN.md` blots
   PKA/AKAP9 as the headline targets with PLK1 as a negative comparator.
3. **Reuse the machinery on a new target.** The validation suite and the AF3 scorecard are
   target-agnostic; point them at another viral bait.
