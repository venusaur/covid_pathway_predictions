# Literature dive: do any open Nsp13 -> PLK1 questions have leads?

A focused sweep of the published literature (2020-2025) against the five open questions. Full
source list with per-source verdicts: `SOURCES.md`. The automated deep-research harness was cut
short by a session rate limit; this is a manual pass and should be re-verified with the harness
after the reset.

## Bottom line

No published study confirms an Nsp13 <-> PLK1 interaction, and a second angle now weakens the
PLK1-specific claim: the field has spent four years validating Nsp13's host partners by co-IP and
imaging, and PLK1 / the centrosome is conspicuously not among them.

- The **upstream premise holds**: SARS-CoV-2 disrupts the host cell cycle, and Nsp13 engages
  cell-cycle and kinase machinery.
- **Nsp13's validated partners are all elsewhere**: STAT1 (IFN suppression), EWSR1 (RNA
  unwinding), the PKA-CREB1 axis (proliferation), and HSP70/90 + ODC (chaperone/polyamine). Each
  confirmed by co-IP and/or immunofluorescence in a primary paper. None is PLK1 or a CEP protein.
- The **PLK1 / centrosome link remains Gordon-AP-MS-only**, never independently reproduced by
  co-IP or imaging, and every well-supported cell-cycle mechanism runs through a *different* axis
  (DNA polymerase delta / S-phase, or the N protein), not PLK1 / G2-M.

One genuinely positive thread: the PKA-CREB1 paper experimentally validates Nsp13-PKA (PRKACA), a
kinase that *was* a Gordon-predicted Nsp13 interactor. So Gordon's Nsp13 predictions can be real
and functional. And PKA is anchored at the centrosome by AKAP9 (itself a Gordon Nsp13 prey), so if
there is any centrosome connection it may run through the *validated* PKA axis rather than PLK1.
That is a sharper, more testable reframing than the original PLK1 guess.

Net: the literature makes the *cell-cycle* premise plausible, the *PLK1-specific* claim less
likely than before, and points at PKA/AKAP9 as the better-supported route to the centrosome if one
exists. The falsifiable core is unchanged: it still needs the co-IP.

## Per-question findings

**Q1 - Any experimental Nsp13 <-> PLK1 interaction?** No, and the negative is now stronger.
No co-IP, AP-MS, BioID, Y2H, or structural study reports it. Several primary papers have
characterized Nsp13's host interactions by co-IP and imaging (STAT1, EWSR1, PKA-CREB1, HSP70/90,
ODC) and none surfaces PLK1 or a centrosome protein. I checked the one lead I had flagged as
high-value, the May et al. BioID proximity proteome (PMC8513853): Nsp13 is not one of its 17
baits, so it does not apply. Stukalov's independent AP-MS recovered only RNF8 for Nsp13. Across
every independent Nsp13-focused study, the PLK1/centrosome axis is absent.

**Q2 - Does the virus disrupt the centrosome / cell cycle / mitosis?** Strong support, wrong axis.
This is well established. SARS-CoV-2 induces S and G2-M arrest (STTT 2022; PMC10972258), the N
protein delays S-phase, and the SARS-CoV Nsp13 ortholog drives S-phase arrest by binding the p125
subunit of DNA polymerase delta and synergizing with Nsp12 (PMC9937006). That last one is the
closest real analogue to our hypothesis: Nsp13 *does* engage a cell-cycle machine. But the
mechanism is DNA-replication / S-phase via Pol delta, not the PLK1-regulated centrosome / G2-M.
So this supports "Nsp13 touches the cell cycle" while pointing away from PLK1 specifically.

**Q3 - Is PLK1 / the centrosome functionally relevant to the virus?** Mild refute.
No SARS-CoV-2 antiviral data exists for PLK1 inhibitors (BI2536, volasertib); their antiviral use
is HIV latency reactivation, via a dual PLK/bromodomain activity (s41598-018-21942-5). This lines
up with our own CRISPR result (PLK1 is not a host-dependency hit) and lowers the prior on PLK1 as
a drug-repurposing target here.

**Q4 - Does Nsp13 localize to the centrosome / bind CEP-family proteins?** Only the original data.
The Nsp13-centrosome links (CEP250, CEP68, AKAP9) all trace back to the Gordon AP-MS and are
interpreted there as microtubule/trafficking hijack, not mitotic. No independent immunofluorescence
places Nsp13 at the centrosome; the one localization study (STTT 2020) puts Nsp13 in the nucleus
and cytoplasm colocalized with splicing speckles (SC35), not centrosomes. No independent
validation of the CEP-scaffold binding that our bridged model rests on.

**Q5 - Is Nsp13 a kinase substrate / does it modulate a host kinase?** Yes, and one is now validated.
Nsp13-PKA (PRKACA) is experimentally confirmed: the PKA-CREB1 axis regulates coronavirus
proliferation through Nsp13 association (PMC11019953), validating a Gordon-predicted Nsp13 kinase
interactor. Nsp13 also binds STAT1 to block its JAK1 phosphorylation (PMC8574307) and engages TBK1
/ TBKBP1 (IFN suppression, Gordon data). So multiple Nsp13-kinase relationships are real. None is
PLK1, and there is no report of Nsp13 as a PLK1 substrate. The validated PKA link matters for our
model: PKA is centrosome-anchored via AKAP9 (a Gordon Nsp13 prey), making PKA/AKAP9 a
better-evidenced route to the centrosome than PLK1.

## What changed for the project

Nothing overturns or confirms the headline, but three concrete updates:

1. The PLK1-specific claim is weaker than before. Four years of Nsp13-focused co-IP/imaging work
   has validated STAT1, EWSR1, PKA-CREB1, HSP70/90, and ODC, and never PLK1 or a CEP protein. A
   real, prominent Nsp13-centrosome function would likely have surfaced by now. It has not.
2. The better-supported route to the centrosome is PKA/AKAP9, not PLK1. Nsp13-PKA is now
   experimentally validated (PMC11019953), and PKA is centrosome-anchored via AKAP9 (a Gordon
   Nsp13 prey). If the co-IP shows any Nsp13-centrosome signal, PKA/AKAP9 is the more likely bridge
   than the PLK1 axis we started with. Worth adding an anti-PKA/AKAP9 blot to the co-IP.
3. The cell-cycle premise holds but points at S-phase too. Nsp13 (SARS-CoV ortholog) drives
   S-phase arrest via DNA polymerase delta, so a real phenotype may be S-phase/DNA-replication as
   much as mitotic. The imaging arm should look at both.

The BioID proximity proteome (PMC8513853) I had flagged as the top computational lead turned out
not to include Nsp13 as a bait, so it does not apply. The co-IP in `EXPERIMENTAL_DESIGN.md` and
the AlphaFold matrix in `ALPHAFOLD_PROTOCOL.md` remain the only ways to settle the PLK1-specific
claim. The literature did not do it for us, and it tilted the odds against PLK1 specifically.
