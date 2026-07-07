# Functional-genomics triangulation — is PLK1 a required SARS-CoV-2 host factor?

**Why this line is different.** Every test up to now measured *interaction* (STRING, IntAct,
HuRI, AP-MS), so all shared the same literature/coverage biases. Genome-wide CRISPR knockout
screens measure *phenotype* — does removing a gene change viral fitness — so they are a fully
**independent** line that escapes that ceiling. Reproducible check: `functional_genomics_check.py`.

**Data (three published genome-wide KO screens, each internally validated by ACE2 ranking #1):**

| Screen | Cell line | Source |
|---|---|---|
| Calu-3 | human lung | Rebendenne/Bruel et al., *Nat. Genet.* 2022 (s41588-022-01110-2) |
| Vero E6 | primate kidney | same paper |
| Bidirectional | (Biering) | Biering et al., *Nat. Genet.* 2022 (s41588-022-01131-x) |

Positive z / positive-selection enrichment = KO **protects** cells = pro-viral host dependency
factor. The residual-z and MAGeCK methods already regress out general gene essentiality.

## Result — PLK1 is not a host factor in any screen

| gene | Calu-3 z (rank / 20,513) | Vero z (rank / 20,928) | Biering pos-FDR (rank / 19,112) |
|---|---|---|---|
| **PLK1** | 1.09 (965) | −0.23 (12,538) | 1.00 (18,382) |
| CEP135 | −0.65 | −0.95 | 0.83 |
| CNTRL | 0.07 | 0.22 | 0.60 |
| PCNT | −0.07 | −0.45 | 1.00 |
| CDK5RAP2 | −0.72 | −0.71 | 1.00 |
| AKAP9 | 0.10 | 0.66 | 0.84 |
| CENPF | 0.50 | −1.63 | 0.93 |
| (validation) ACE2 | z 8.05, **#1** | z 11.7, **#1** | FDR 3e-4, **#1** |

Hit thresholds in these screens: |z| ≥ 2 (only 92 genes in Calu-3, 368 in Vero) or FDR < 0.1.
**Neither PLK1 nor any centrosome-module gene clears them, in any screen.** PLK1's best showing
is Calu-3 z = 1.09 (top ~5%, but sub-threshold); it is neutral in Vero and non-significant in
Biering. The module scatters around zero.

## Interpretation — a negative for *requirement*, not for *interaction*

Three points, in order of importance:

1. **This does NOT refute the Nsp13–PLK1 interaction.** Physical association and fitness
   requirement are different questions. A protein the virus binds or perturbs need not be
   *essential for replication* — it may be redundant, buffered in these cell lines, or relevant
   to pathology rather than to the culture-replication readout the screens select on. Binding ≠
   requirement; a CRISPR-null result constrains function, it does not erase a contact.

2. **It reframes the biological hypothesis productively.** Combined with the motif analysis (no
   direct PLK1 docking motif → bridged model), the cleanest reading is: *if* Nsp13 engages the
   PLK1-regulated centrosome, the consequence is more likely **cell-cycle / centrosome
   perturbation (a pathogenesis phenotype)** than a replication-essential host-factor role. That
   is a sharper, differently-testable claim (look at centrosome number, mitotic arrest, γ-tubulin
   recruitment in infected cells — already in the wet-lab plan) than "PLK1 is needed to replicate."

3. **It tempers the drug-repurposing angle.** The earlier note floated PLK1 inhibitors
   (BI2536, volasertib) as candidate antivirals. Three screens saying PLK1 is not a host
   dependency factor lowers that prior substantially — and PLK1 inhibitors are anti-mitotic
   cytotoxics, so any apparent antiviral signal would be badly confounded by toxicity. Not a
   priority lead.

## The honest confound, handled

PLK1 is essential for mitosis, so naively its sgRNAs should deplete (toxicity) and look like a
false "anti-viral" hit. Two things rule that out here: the screens use essentiality-corrected
scores (residual-z, MAGeCK), and PLK1 came out **neutral-to-slightly-positive**, not depleted —
so the null is a real null, not toxicity masking a signal.

## Where this leaves the hypothesis

Convergent with everything else: Nsp13–PLK1 is a specific, cross-method-reproducible *interaction*
hypothesis (test 12) that is **not** a functional host-dependency relationship (this analysis).
The most likely biology, if real, is centrosome/cell-cycle disruption during infection — which is
exactly what the co-IP + imaging experiment in `EXPERIMENTAL_DESIGN.md` is built to detect, and
what the AlphaFold bridged-vs-direct matrix in `ALPHAFOLD_PROTOCOL.md` is built to localize.
