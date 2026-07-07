# Experimental design note — testing SARS-CoV-2 Nsp13 ↔ human PLK1

**Purpose.** Every computational avenue has now been exhausted (14 validation tests; see
`validation/VALIDATION_REPORT.md`). Nsp13→PLK1 is a specific, cross-method-reproducible
hypothesis (STRING + IntAct physical graph, p=0.0005) that **no available dataset can confirm
bias-free** — the one unbiased binary interactome (HuRI) is blind to the centrosome. Resolving
it now requires a wet-lab experiment. This note specifies that experiment.

**The precise hypothesis, and the key ambiguity to resolve.** The evidence is that PLK1 is
tightly co-functional with the centrosome/centriole proteins Nsp13 directly binds
(CEP135, CEP43, CNTRL, NIN, NINL — all confirmed Nsp13 preys and all IntAct physical partners
of PLK1). Sequence analysis (`computational_followups/motif_scan`) finds **no canonical PLK1
Polo-box docking motif (S-[pS/pT]-P) in Nsp13** and only background-level kinase-consensus
sites. Together this favours a **bridged/shared-complex model** — Nsp13 associates with the
PLK1-regulated centrosome module, plausibly *through* the CEP scaffolds rather than by direct
Nsp13–PLK1 contact. The experiment must therefore test not just *whether* they associate but
*whether the association is direct or bridged*.

---

## Primary experiment — reciprocal co-immunoprecipitation

**System.** Use the Gordon et al. platform for comparability: 2×Strep-tagged Nsp13 expressed in
HEK-293T. Then confirm in a physiologically relevant, infectable lung line (A549-ACE2 or
Calu-3), and ideally under **authentic SARS-CoV-2 infection** with an anti-Nsp13 antibody, to
rule out overexpression artifacts.

**Reciprocal pulldowns (both directions — reciprocity is the single biggest false-positive
filter):**
1. **Forward:** Strep-Nsp13 affinity pulldown → immunoblot for **endogenous PLK1**.
2. **Reverse:** endogenous PLK1 IP (or GFP-PLK1) → immunoblot for Nsp13.

**Internal positive controls (ties the assay to the computational evidence).** On the same
Nsp13 pulldown, blot for **CEP135, CNTRL, NINL, CEP250** — the MiST-confirmed Nsp13 preys that
are also PLK1 partners. These *must* co-purify; if they don't, the pulldown failed and a
negative PLK1 result is uninterpretable.

**Negative controls.** (a) Strep-only / empty vector; (b) an unrelated tagged viral bait that
should *not* pull PLK1 — **Nsp8 or ORF9b**; (c) IgG isotype control for the reverse IP.

**Cell-cycle staging — do not skip.** PLK1 abundance and activity peak in G2/M. Run the co-IP in
**both asynchronous and G2/M-synchronized** cells (thymidine–nocodazole release). A cell-cycle-
gated interaction will be invisible in asynchronous lysate — a plausible false-negative mode.

**Phospho-dependence arm.** Although Nsp13 lacks a canonical PBD motif, include a
**λ-phosphatase-treated lysate** arm. Loss of the interaction on dephosphorylation would
indicate a phospho-primed (PBD-type) contact after all; insensitivity is consistent with the
bridged model.

---

## Disambiguating **direct vs bridged** (the crux)

Cellular co-IP cannot tell a direct Nsp13–PLK1 contact from bridging via a centrosome scaffold.
Two orthogonal tests settle it:

- **In vitro reconstitution.** Purified recombinant Nsp13 + PLK1 (His/GST pulldown; quantify
  affinity by BLI or SPR). **Binding here = direct.** No binding, despite a positive cellular
  co-IP, points to bridging.
- **Bridge-dependence in cells.** Repeat the cellular co-IP after **siRNA knockdown of the
  candidate bridge (CEP135 / CNTRL)**. If the Nsp13–PLK1 co-IP collapses, the association is
  scaffold-mediated.

| Forward+reverse co-IP | Recombinant binding | Bridge-KD sensitive | Interpretation |
|---|---|---|---|
| + | + | (n/a) | **Direct Nsp13–PLK1 interaction** |
| + | − | yes | **Bridged** through the centrosome scaffold |
| + | − | no | Indirect via an untested third partner |
| − | − | — | Not supported (note co-IP false-negative modes: transient / cell-cycle-gated) |

---

## Orthogonal validation

- **Proximity labeling.** TurboID- or BioID-Nsp13 → streptavidin capture → MS. Expect PLK1 and
  the centrosome module enriched in the Nsp13 proximity interactome.
- **Imaging.** Immunofluorescence — does Nsp13 localize to the centrosome (co-stain γ-tubulin /
  pericentrin, PLK1)? Live or infected cells.

## Functional follow-up (only if association confirmed)

- **Antiviral relevance (ties back to the drug-repurposing theme of the source paper).** PLK1
  inhibitors — **BI2536, volasertib, or onvansertib** (clinical-stage) — dose–response vs
  SARS-CoV-2 replication (plaque assay / RT-qPCR) in Calu-3 or A549-ACE2, with a cytotoxicity
  counter-screen. PLK1 siRNA/CRISPR as genetic confirmation.
- **Mechanism.** Does Nsp13 perturb centrosome function — centrosome number, γ-tubulin
  recruitment, mitotic index? Is Nsp13 a PLK1 *substrate* (in-vitro kinase assay ± PLK1) or a
  *regulator/hijacker*? The weak motif evidence predicts regulator, not substrate.

## Rigor

Biological triplicate; pre-registered analysis; blinded quantification for imaging; report
negatives. The internal centrosome-prey positive controls are the linchpin — they convert a
bare "we saw a band / we didn't" into an interpretable result.

---

*Prepared alongside the computational validation in this repository. Ready-to-run structural-
prediction input and the sequence-motif analysis that motivated the direct-vs-bridged framing
are in `experimental/` and `experimental/COMPUTATIONAL_PATHWAYS.md`.*
