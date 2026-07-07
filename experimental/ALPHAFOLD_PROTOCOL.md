# AlphaFold-Multimer protocol — structural test of Nsp13 ↔ PLK1

**Why this is the right test.** Every interaction-database avenue hit the shared-study-bias
ceiling (validation tests 8–14). AlphaFold-Multimer predicts physical docking from **sequence +
structure alone, with zero dependence on any interaction database** — so it is the one remaining
computational line that is both bias-free and able to give *positive* evidence.

**Honest status.** This was **not run in the analysis environment** — AlphaFold-Multimer needs a
GPU and hundreds of GB of MSA databases. What is prepared here is everything needed to run it
rigorously elsewhere (ColabFold or the AF3 server), plus the pre-registered interpretation
framework so the result can't be talked into meaning whatever we want. Inputs:
`experimental/alphafold_inputs/*.fasta`; manifest: `run_matrix.csv`.

---

## The run matrix (7 jobs) — designed, not a single naive run

Running only "Nsp13 + PLK1 full-length" would be uninterpretable: two large multidomain proteins
can find a spurious interface, and there is no baseline to judge the score against. The matrix
fixes both — it **localizes** the interface and it **calibrates** the score with controls.

| Run | Complex | Role | Question it answers |
|---|---|---|---|
| 01 | Nsp13 + PLK1 (full) | test | Can they dock at all? |
| 02 | Nsp13 + PLK1-**PBD** (367–603) | test | Via PLK1's canonical **partner-docking** Polo-box domain? |
| 03 | Nsp13 + PLK1-**kinase** (1–345) | test | Or the kinase domain (substrate-like)? |
| 04 | Nsp13 + **CEP135** | **positive / bridge** | Nsp13's *confirmed* direct prey — AF should find this interface |
| 05 | Nsp10 + Nsp16 | **positive / pipeline** | Known SARS-CoV-2 heterodimer (solved) — proves the run finds true interfaces |
| 06 | Nsp13 + Catalase | **negative / null** | Unrelated protein — sets the false-positive floor |
| 07 | Nsp13 + Pyruvate kinase | **negative / null** | Second decoy — the structural null baseline |

Domain boundaries are real (UniProt): PLK1 kinase 53–305, Polo-box domain 410–592 (phospho-
peptide pocket at 538–540); Nsp13 zinc-binding domain 1–113, helicase ATP-binding 257–438.

---

## How to run

**Option A — ColabFold (free GPU, AlphaFold2-Multimer).** For each FASTA: open the ColabFold
"AlphaFold2_mmseqs2" notebook, paste the two-chain sequence (chains separated by `:` or as the
multi-record FASTA), set `model_type = alphafold2_multimer_v3`, `num_recycles = 6` (bump higher
for the large CEP135 job), enable `--amber` relaxation. Download the ZIP.

**Option B — AlphaFold3 server** (alphafoldserver.com): paste the two sequences as separate
protein entities. AF3 also lets you add a **phosphate PTM** — relevant here (see caveat 2).

Large jobs (04, 06, 07 exceed 1100 residues) run fine on Colab Pro or a local A100.

---

## What to record (per run, best model)

- **ipTM** — interface predicted TM-score. The primary readout.
- **pTM** — global.
- **Interface PAE** — is there a *compact, low-PAE block* between the two chains? This matters
  more than ipTM alone: a low-PAE contiguous interface is real signal; a high-PAE contact is not.
- **Interface pLDDT** and the **specific contact residues** (for the positive case → mutagenesis).

`experimental/alphafold_inputs/run_matrix.csv` has a row per job to fill in.

---

## Pre-registered interpretation (decide the thresholds *before* looking)

**Step 1 — validate the pipeline with the controls, or trust nothing:**
- Run 05 (Nsp10+Nsp16) **must** score high (ipTM ≳ 0.8, low-PAE interface). If it doesn't, the
  MSA/run is broken — stop and fix before interpreting anything.
- Runs 06 & 07 (decoys) define the **null band** — the ipTM non-binders score on this hardware
  (typically ≲ 0.3). This is the AlphaFold analog of the degree-preserving null in test 6.

**Step 2 — read the tests against that calibrated baseline:**

| Pattern across runs 01–04 | Interpretation | Feeds wet-lab (`EXPERIMENTAL_DESIGN.md`) |
|---|---|---|
| PBD (02) high & low-PAE, esp. near the 538–540 pocket | **Direct binding via the canonical PBD** | Prioritize in-vitro reconstitution; mutate the interface residues |
| Kinase (03) high, PBD low | Direct, substrate-like contact | In-vitro kinase assay + binding |
| All PLK1 runs (01–03) at/below the decoy null, **but CEP135 (04) high** | **Bridged model** — no direct interface; association is scaffold-mediated | Prioritize the CEP135-knockdown arm; expect scaffold-dependent co-IP |
| Full (01) high but both domains (02,03) at null | Avidity/artifact — interpret with caution | Re-run with more recycles/seeds before trusting |
| CEP135 (04) *also* at null | Pipeline can't model Nsp13 interfaces at all | AF is uninformative for this bait; rely on the wet-lab |

This is concordant-by-design with the motif result: Nsp13 has **no** canonical PLK1 Polo-box
docking motif, which already predicts the **bridged** row. A confident positive PBD interface
would be the surprising, high-value outcome — it would overturn the motif-based expectation and
hand the experiment exact residues to test.

---

## Caveats — stated up front, not after a convenient result

1. **Host–virus reliability.** AF-Multimer is calibrated mainly on co-evolved intra-species
   complexes; virus–host pairs have shallow paired MSAs, so confidence is systematically lower
   and a *negative* result does **not** refute the hypothesis.
2. **Phospho-priming.** PLK1's PBD binds *phosphorylated* partners. Standard AF models the
   unmodified sequence and can therefore **miss** a phospho-dependent interface — another reason
   a negative PBD run is inconclusive. (AF3's phosphate PTM partially addresses this; worth a
   variant run.)
3. **A positive is much stronger than a negative here.** Given 1 and 2, treat a high-confidence,
   low-PAE interface as strong bias-free support; treat a null as "uninformative," not "refuted."
4. **Validate against experiment.** Nsp13 has solved structures (PDB 6ZSL, 7NIO, 7NNG) and PLK1's
   PBD is solved (3D5U–3DB8) — compare AF monomer domains to these before trusting the complex.

---

### Bottom line

The matrix turns one ambiguous prediction into a calibrated experiment with its own null and
positive controls — the same standard the rest of this project held. It cannot *confirm*
Nsp13–PLK1 on its own (a wet-lab co-IP remains decisive), but a confident PBD interface would be
the strongest bias-free evidence obtainable computationally, and either way it sharpens the
direct-vs-bridged question the bench experiment is built to resolve.
