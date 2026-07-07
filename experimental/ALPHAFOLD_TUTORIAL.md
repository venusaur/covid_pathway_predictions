# AlphaFold tutorial — step by step, what to actually do

A hands-on walkthrough for running the Nsp13↔PLK1 structural test, written to be followable with
no prior AlphaFold experience. (The *scientific* design — why these 7 jobs, how to interpret the
scores — is in `ALPHAFOLD_PROTOCOL.md`. This file is the button-pushing.)

**Time & cost:** ~15–60 min per job on a free Colab GPU; the whole 7-job matrix in an afternoon.
Free. No install.

**What you're doing in one sentence:** giving AlphaFold two protein sequences and asking "if these
two proteins were in a tube together, would they stick, and where?"

---

## Before you start — the inputs are already made

Everything you paste is in **`experimental/alphafold_inputs/`** (7 files) — you don't write any
sequences yourself. The jobs, in the order to run them:

| File | What it tests | Run it… |
|---|---|---|
| `05_nsp10__nsp16_POSCTRL_pipeline.fasta` | **Run this FIRST** — a known real complex; proves the tool works | 1st |
| `06_nsp13__CAT_NEGCTRL.fasta` | a non-binder — sets your "looks like nothing" baseline | 2nd |
| `07_nsp13__PKM_NEGCTRL.fasta` | second non-binder | 3rd |
| `01_nsp13__plk1_full.fasta` | the actual question | 4th |
| `02_nsp13__plk1_PBD.fasta` | is the contact via PLK1's docking domain? | 5th |
| `03_nsp13__plk1_kinase.fasta` | …or its kinase domain? | 6th |
| `04_nsp13__cep135_POSCTRL_bridge.fasta` | Nsp13's known partner — the "bridged?" test | 7th |

> **Why run the controls first?** If the known complex (job 05) doesn't score high, something is
> wrong with the run and *none* of your other numbers mean anything. And the two non-binders
> (06, 07) tell you what a "no" looks like *on your hardware*, so you can judge PLK1 against a real
> baseline instead of a number you guessed. This is the whole reason the matrix exists.

---

## Route A — ColabFold (recommended, free)

1. **Open the notebook.** Go to **https://github.com/sokrypton/ColabFold** and click the
   **"AlphaFold2_mmseqs2"** Colab badge (opens Google Colab in your browser). You need a Google
   account; nothing to install.
2. **Turn on the GPU.** Colab menu → *Runtime* → *Change runtime type* → Hardware accelerator =
   **T4 GPU** → Save.
3. **Enter the sequence.** In the `query_sequence` field, paste the **two sequences joined by a
   colon `:`** — the colon tells AlphaFold "these are two separate chains." Open a FASTA file from
   `alphafold_inputs/`, and paste the two sequence blocks separated by `:` like:
   ```
   AVGACVLCNSQ...(Nsp13)...ATLQ:MSAAVTAGKLARAP...(PLK1)...
   ```
   Set `jobname` to something you'll recognize, e.g. `05_nsp10_nsp16`.
4. **Settings that matter** (leave the rest default):
   - `num_recycles` = **6** (for the big jobs 04/06/07, use 12 if it's stable)
   - `template_mode` = none is fine
   - tick **`use_amber`** (cleans up the final model)
5. **Run it.** *Runtime* → *Run all*. It will (a) build an MSA via the MMseqs2 server, then
   (b) predict 5 models. Watch for the green checkmarks.
6. **Download the results.** When done it auto-downloads (or offers) a **`.zip`**. Inside:
   - `*_scores_rank_001_*.json` — **the numbers you need** (`iptm`, `ptm`, `pae`)
   - `*_predicted_aligned_error_v1.json` and the PAE PNG — the interface-error map
   - `*_relaxed_rank_001_*.pdb` — the 3-D model (open in Mol* / PyMOL / ChimeraX)
7. **Repeat** for all 7 files.

## Route B — AlphaFold3 server (no Colab, and can model phosphorylation)

1. Go to **https://alphafoldserver.com** and sign in.
2. Add **two "Protein" entities**; paste the Nsp13 sequence into one and PLK1 into the other
   (from `01_nsp13__plk1_full.fasta`). No colon here — they're separate boxes.
3. *(Optional but valuable)* AF3 lets you add a **phosphorylation PTM**. PLK1's docking domain
   only grabs *phosphorylated* partners, so a standard run can miss a real contact — see caveat 2
   below. A phospho-variant run is a nice extra.
4. Submit; read `ipTM`/`pTM` and the PAE plot in the results view.

---

## Reading the result — a 30-second version

For each job, open the scores JSON and note three things:

- **ipTM** (0–1): interface confidence. **This is the headline number.**
- **pTM** (0–1): whole-model confidence (less important here).
- **Interface PAE**: in the PAE plot, look at the **off-diagonal block** (chain-A rows vs chain-B
  columns). **Dark/low there = a confident interface.** Light/yellow = the two chains are just
  floating near each other — not real.

Then judge against your own controls, not a textbook cutoff:

```
Did job 05 (Nsp10+Nsp16) score high (ipTM ~0.8+, dark off-diagonal)?
   NO  -> stop. The run is broken; fix before trusting anything.
   YES -> continue.

What ipTM did the decoys (06, 07) give?  Call that your NULL (often ~0.2-0.3).

Now look at PLK1 (jobs 01/02/03):
   ipTM clearly ABOVE null AND a dark off-diagonal PAE block
        -> real predicted interface. Note WHICH PLK1 domain (02 vs 03) and which residues.
   ipTM at/below the null, PAE all light
        -> no predicted direct interface.

Cross-check with job 04 (Nsp13+CEP135):
   CEP135 high but PLK1 low  -> supports the BRIDGED model (matches the motif result).
   Both low                  -> AF can't model Nsp13 interfaces; rely on the wet-lab.
```

Fill your numbers into `experimental/alphafold_inputs/run_matrix.csv` as you go.

---

## Three caveats — read them before over-reading a result

1. **A "no" is weak; a "yes" is strong.** AlphaFold is trained mostly on proteins from the *same*
   organism that evolved together. Virus–host pairs have thin evolutionary signal, so a low score
   does **not** disprove the interaction — it's just uninformative. A *high*, low-PAE interface,
   though, is strong bias-free evidence.
2. **PLK1's docking domain needs a phosphate.** PLK1's Polo-box domain only binds partners that
   are already phosphorylated. A standard AlphaFold run uses the plain sequence and can therefore
   *miss* that contact. If job 02 comes back low, that's expected, not conclusive — try the AF3
   phospho-variant (Route B, step 3).
3. **Check the monomers look right first.** Nsp13 and PLK1 have known experimental structures
   (Nsp13: PDB **6ZSL**; PLK1 Polo-box: PDB **3D5U**). If AlphaFold's single-chain fold for either
   protein looks nothing like these, don't trust the complex.

---

## What a good outcome looks like (worked example)

Say you get: job 05 ipTM 0.86 (pipeline works ✓), jobs 06/07 ipTM 0.24 and 0.21 (null ≈ 0.22),
job 01 ipTM 0.31, job 02 (PBD) ipTM 0.29, job 03 (kinase) ipTM 0.27, job 04 (CEP135) ipTM 0.71
with a dark PAE block.

**Read:** PLK1 runs are at the null → no *direct* Nsp13–PLK1 interface predicted; but Nsp13 docks
its known partner CEP135 well → AlphaFold *can* find Nsp13 interfaces, so the PLK1 null is
meaningful, not a pipeline failure. **Conclusion:** supports the **bridged** model — take it to the
bench with the CEP-knockdown arm (see `EXPERIMENTAL_DESIGN.md`). *(Illustrative numbers — you'll
generate the real ones.)*

The opposite outcome — a confident, low-PAE PLK1-PBD interface — would be the surprising,
high-value result: direct binding, with exact residues to mutate. Either way, the run is
decision-shaped, not a coin toss.
