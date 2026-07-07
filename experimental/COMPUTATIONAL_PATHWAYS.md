# Remaining computational pathways for Nsp13 → PLK1

The 14-test validation suite exhausted the **interaction-database** avenues: every dataset that
covers the centrosome is literature-curated, and the one unbiased binary map (HuRI) can't see it
(test 14). So *association-based* computation has hit its ceiling. But several **structurally
different** computational lines remain — ones that do not query interaction databases at all and
so are immune to the shared-study-bias problem. Ranked by strength and independence:

| # | Pathway | Independent of PPI-database bias? | Status here |
|---|---|---|---|
| 1 | **AlphaFold-Multimer / AF3 structural docking** | **Yes — biophysical, sequence-only** | **7-job run matrix + protocol prepared** (`ALPHAFOLD_PROTOCOL.md`) |
| 2 | Ternary / bridged-topology modelling | Yes | Specified |
| 3 | PLK1 sequence-motif analysis | Yes | **Done** (below) |
| 4 | Functional-genomics triangulation (CRISPR / drug screens) | Yes — phenotype, not interaction | **Done** — 3 screens; PLK1 not a host factor (`FUNCTIONAL_GENOMICS.md`) |
| 5 | Cross-coronavirus interactome conservation | Partly | Specified (existing data too sparse) |

---

## 1. AlphaFold-Multimer / AF3 — the strongest remaining test  ·  *input ready*

This is the one that genuinely escapes the ceiling: it predicts whether Nsp13 and PLK1 can
**physically dock** from sequence + structure alone, with **zero dependence on any interaction
database**. No literature bias, no Y2H coverage gap.

- **Fully prepared:** a **7-job run matrix** — not one naive full-length run — with domain-split
  tests (PLK1 kinase vs Polo-box domain), a positive/bridge control (Nsp13 + its confirmed prey
  CEP135), a pipeline positive control (Nsp10 + Nsp16), and two structural-null decoys to
  calibrate the ipTM baseline. Inputs in `experimental/alphafold_inputs/`; full method,
  metrics, and a **pre-registered interpretation decision-tree** in `ALPHAFOLD_PROTOCOL.md`.
- **How to run:** ColabFold (AlphaFold2-Multimer) or the AF3 server — both take the two-chain
  FASTAs directly. No local GPU needed via Colab; not runnable in the analysis sandbox.
- **Read the right metrics:** interface **ipTM** judged *against the decoy null band*, and the
  **interface PAE** (a compact low-PAE block is the real signal). Validate with the controls
  before trusting any test row.
- **Caveats (up front):** AF-Multimer is less reliable for host–virus pairs, and PLK1's PBD binds
  *phospho*-primed partners that a standard run won't model — so a **negative does not refute**;
  a confident low-PAE interface would be strong bias-free support and would name the residues to
  mutate. The motif result already predicts the *bridged* outcome, making a positive PBD interface
  the surprising, high-value result to watch for.

## 2. Ternary / bridged-topology modelling — tests the direct-vs-bridged question *in silico*

The motif analysis (§3) predicts the association is **bridged** through a CEP scaffold rather than
a direct Nsp13–PLK1 contact. Structure prediction can pre-test this before the bench does:

- Model **Nsp13 + CEP135** and **CEP135 + PLK1** pairwise, and the **Nsp13 + CEP135 + PLK1**
  ternary. If the two scaffold interfaces score well but the direct Nsp13–PLK1 interface does
  not, that *computationally* favours bridging — and tells the wet-lab experiment (see
  `EXPERIMENTAL_DESIGN.md`) to expect a scaffold-dependent, not direct, co-IP.
- Cheap, fully in-scope, and directly shapes the experimental interpretation table.

## 3. PLK1 sequence-motif analysis  ·  **done** (`experimental/motif_scan.py`)

PLK1 recognises partners two ways; we scanned Nsp13 for both:

- **Polo-box docking motif** `S-[pS/pT]-P` (Elia et al. 2003) — the canonical stable-recruitment
  motif: **0 found** in Nsp13 (chance expectation ≈0.3, so absence is only mildly informative).
- **Kinase-substrate consensus** `[D/E/N]-X-[S/T]-Φ` — **5 sites** (positions 142, 197, 349,
  375, 383), i.e. **background rate**, no enrichment.

**Interpretation:** no evidence for a direct kinase–substrate or PBD-docking relationship →
favours the **bridged/shared-complex model**. This reframes the hypothesis productively: Nsp13
engages the **PLK1-regulated centrosome module** (via the CEP/centriolar scaffolds it directly
binds), rather than PLK1 itself — a sharper, more testable claim than "Nsp13 binds PLK1."

## 4. Functional-genomics triangulation — **done**  (`functional_genomics_check.py`, `FUNCTIONAL_GENOMICS.md`)

Checked PLK1 and the centrosome module across **three genome-wide CRISPR KO screens** (Calu-3
lung, Vero E6, and the Biering bidirectional screen — each validated by ACE2 ranking #1).
**PLK1 is not a significant host factor in any of them** (best case Calu-3 z=1.09, top ~5% but
sub-threshold; neutral in Vero; non-significant in Biering), and neither is any centrosome gene.

Crucially, this is a negative for *functional requirement*, **not** for the interaction: binding
≠ fitness-essentiality. It productively reframes the hypothesis — if Nsp13 engages the
PLK1-centrosome module, the likely consequence is **cell-cycle/centrosome perturbation
(pathogenesis)** rather than an essential replication-factor role — and it tempers the
PLK1-inhibitor drug-repurposing angle (also confounded by anti-mitotic cytotoxicity). The
essentiality confound is handled: the screens use essentiality-corrected scores and PLK1 came
out neutral, not depleted.

## 5. Cross-coronavirus interactome conservation — currently under-powered

If Nsp13's centrosome hijacking is conserved in SARS-CoV / MERS, that's evolutionary support.
Note the honest limitation: the one comparative dataset to hand (Stukalov et al. 2021, which
profiled SARS-CoV Nsp13) is **too sparse for Nsp13** — its SARS-CoV-2 Nsp13 recovered a single
interactor and its SARS-CoV Nsp13 hit list shows no centrosome cluster. A real test needs a
dedicated, well-powered SARS-CoV / MERS Nsp13 AP-MS, not reanalysis of existing sparse data.

---

### Bottom line

Two of these pathways are now executed (motif scan §3, functional genomics §4) and they
converge on the same refined picture: Nsp13–PLK1 is a specific *interaction* hypothesis, most
likely **bridged through the CEP centrosome scaffolds** (no direct docking motif) and
**not a functional host-dependency relationship** (not a CRISPR hit in three screens) — so its
biology, if real, is centrosome/cell-cycle perturbation, not a replication-essential factor.

The **AlphaFold-Multimer run (7-job matrix prepared) is the single highest-value remaining
computational step** — the only avenue that is both bias-free and able to give *positive*
structural evidence, and it is now designed to directly test the bridged-vs-direct question the
other two pathways raised. None of it substitutes for the wet-lab co-IP in
`EXPERIMENTAL_DESIGN.md`, but together they have turned a vague "Nsp13 binds PLK1" into a
precise, mechanism-specific, experimentally-actionable model.
