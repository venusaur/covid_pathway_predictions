# Validation report: does the missing-component prediction actually work?

This report runs 14 adversarial tests against the two headline claims from the earlier analysis:
**Nsp13 -> PLK1** and **Nsp4 -> TOMM70** (both "predicted missing components" via STRING
guilt-by-association). Tests 1-8 probe the two specific claims; tests 9-14 step back and ask
whether the prediction *method itself* has any real predictive value against independent data
and independent network substrates (it largely does not -- see tests 10-11 -- with the notable
exception that Nsp13->PLK1 reproduces on a physically-constructed graph, test 12; test 14 then
maps the ceiling on how far computation can confirm it at all). Every test below is implemented
as a runnable script in this folder and was executed against live data -- nothing here is
asserted without a number behind it.

## Summary table

| # | Test | Script | Verdict |
|---|---|---|---|
| 7 | Identifier integrity | `01_identifier_integrity.py` | **PASS** -- 331/332 exact UniProt/HGNC matches, 1 legitimate synonym, zero broken joins |
| 4 | Coverage/denominator audit | `02_coverage_audit.py` | **CONCERN CONFIRMED** -- Reactome covers only 68.7%, KEGG 46.4% of our 332 genes; coverage is bait-dependent |
| 5 | Strip circular annotations | `03_strip_circular.py` | **PASS** -- removing 5 SARS-named Reactome terms doesn't reorder the remaining ranking |
| 2 | Negative controls (E, M, Orf7a) | run inline, see below | **PASS** -- PLK1/TOMM70 never leak into unrelated baits; predictions stay topically anchored |
| 6 | Threshold robustness | `04_threshold_sweep.py` | **SPLIT** -- Nsp13->PLK1 robust across most of the parameter band; Nsp4->TOMM70 fragile, survives only at the loosest confidence setting |
| 1 | Degree-preserving null model | `05_null_model.py` | **PASS** -- 0/25 hits for both PLK1 and TOMM70 under randomized same-size gene sets; not generic hubs |
| 3 | Precision/recall via hold-out | `06_holdout_precision_recall.py` | **REAL BUT MODEST** -- 11.1% recall / 3.5% precision (real) vs 0%/0% (null), but highly bait-dependent (0% for M, Nsp7, Nsp12, Orf9b) |
| 8 | External validation (IntAct) | `07_external_validation.py` | **NOT CONFIRMED** -- TOMM70 is well-established as a SARS-CoV-2 target across 7 independent papers, but always with ORF9b, never Nsp4. PLK1 has 1 independent SARS-CoV-2 record, with Spike, never Nsp13 |
| 9 | Independent AP-MS dataset crosscheck (Stukalov et al. 2021) | `08_independent_dataset_crosscheck.py` | **NOT SUPPORTED, LOW POWER** -- Stukalov's own published hit list (not a third-party curation) reports exactly 1 significant Nsp13 interactor (RNF8), no PLK1. But Nsp13 is the single sparsest bait in their entire SARS-CoV-2 panel (dead last of 22 baits) -- their assay recovered almost nothing for this bait, so this is a weak negative, not a strong one |
| 10 | Precision/recall against genuinely independent ground truth | `09_independent_ground_truth_precision_recall.py` | **METHOD MOSTLY FAILS TO GENERALIZE** -- using Gordon et al.'s full real prey lists (no holdout) to predict what Stukalov et al. independently found, across 18 testable baits: mean recall 0.07% (real) vs 0.02% (degree-matched null) -- both near zero. The one exception, Orf3a, recovered 4/329 novel Stukalov interactors (VPS16, VPS18, VPS41, VPS45 -- all HOPS complex subunits, matching real, independently-published ORF3a-HOPS/autophagy biology), against 1/329 for its own null. This is a materially harder and more honest test than test 3's internal holdout, and it substantially lowers confidence that this guilt-by-association method generalizes to real, undiscovered biology beyond the training dataset it was tuned on |
| 11 | Recall-gap diagnosis (why some baits work) | `10_recall_gap_diagnosis.py` | **DIAGNOSIS -- SUCCESS IS LITERATURE-BIASED** -- internal-holdout recall tracks how densely a bait's known preys are already interconnected in STRING (Spearman rho +0.65 for edge density), confirming the method only works where the bait's biology is already a tight cluster. But the single strongest predictor of recall is the *text-mining* channel (rho +0.78), not hard experimental evidence (escore rho +0.50, below the n=10 significance threshold). Nsp8 is the clean counterexample: richest in experimental evidence yet ~0 recall, because its interactors aren't the literature-co-mentioned ones. This mechanistically explains test 10: the method's apparent successes ride on literature co-mention bias, which by construction does not transfer to novel, less-studied interactions |
| 12 | Alternative-substrate reproduction (IntAct physical graph) | `11_alternative_substrate_intact.py` | **SPLIT -- PLK1 REPRODUCES, TOMM70 DIES** -- redoing guilt-by-association on IntAct's curated *physical* interactions (no text-mining channel, a structurally different graph) instead of STRING. **Nsp13->PLK1 reproduces**: 6 of 40 Nsp13 preys (CEP135, CEP43, CNTRL, NIN, NINL, TBK1 -- all centrosome/centriole) are direct physical partners of PLK1, vs null mean 1.3, empirical p=0.0005; PLK1 ranks in the top ~1.2% of IntAct candidates, surrounded by other centrosome proteins. **Nsp4->TOMM70 does not**: 0 of 8 preys overlap TOMM70's IntAct partners, unranked. This is the first positive cross-method evidence for PLK1 (it is not a STRING/text-mining artifact) and the final nail for TOMM70. Caveat: IntAct is also literature-curated, so it shares *some* study bias, but the graph is built from physical binding, not co-mention |
| 13 | Broadened pathway coverage | `12_broaden_coverage.py` | **GAP IS STRUCTURAL** -- adding WikiPathways 2023, BioPlanet 2019, and MSigDB Hallmark 2020 lifts union coverage of the 332 genes only from 95.2% to 96.4%; 12 of the 16 genes unannotated in GO/KEGG/Reactome stay dark across all six independent libraries, and the same membrane/small-ORF baits (Orf7a, Spike, Nsp12, E) stay coverage-poor. The enrichment blind spots reflect genuinely under-characterized biology, not one database's idiosyncrasy -- so cross-bait enrichment comparisons remain not-on-equal-footing no matter which library you pick |
| 14 | Non-curated binary interactome (HuRI) | `13_binary_interactome_huri.py` | **UNINFORMATIVE -- COVERAGE CEILING, not a refutation** -- the reviewer asked for an unbiased dataset that could actually confirm PLK1, since STRING/IntAct/Gordon are all literature-curated. HuRI (Luck et al. 2020, systematic all-by-all Y2H, not literature-assembled) is that dataset -- but it is structurally blind to the biology in question: **10 of 14 core centrosome proteins are entirely absent from HuRI**, PLK1 itself has degree 1, and only 1 of the 6 IntAct preys that support PLK1 even appears. A 0-overlap here is a Y2H coverage failure (large coiled-coil scaffolds express poorly in two-hybrid), not evidence against the hypothesis. **The honest ceiling: every dataset that covers the centrosome is literature-curated, and the one unbiased binary map cannot see it -- so the shared-study-bias caveat on Nsp13->PLK1 is not removable with existing computational data. Only a targeted wet-lab co-IP can close it** |

## What this means for the two headline claims

**Nsp13 -> PLK1**: the strongest-standing of the two, and the evidence is genuinely mixed rather
than uniformly negative. In favour: it survives every internal test (null model, negative
controls, most of the threshold band); Nsp13's holdout recall is one of the best in the dataset
(27.5%); and -- most importantly -- it **reproduces on a structurally independent network**
(test 12). Rebuilding the prediction on IntAct's curated *physical*-interaction graph, which has
no text-mining channel, still places PLK1 among Nsp13's centrosome preys: 6 of 40 preys are
direct physical partners of PLK1 (p=0.0005), and PLK1 sits in the top ~1.2% of candidates
surrounded by other centrosome proteins. So the call is not a STRING/text-mining artifact.
Against: no independent SARS-CoV-2 *viral* interactome has directly observed PLK1 binding Nsp13
-- IntAct's cross-database search (test 8) links PLK1 to Spike, not Nsp13, and Stukalov et al.'s
own AP-MS (test 9) reports just one Nsp13 interactor (RNF8, not PLK1), though that miss is
low-power since their Nsp13 recovered almost nothing at all. The honest read: *a specific,
non-artifactual, cross-method-reproducible hypothesis about which human machinery Nsp13 engages
(the PLK1-regulated centrosome), not yet directly confirmed as a physical Nsp13-PLK1 contact.*
A co-IP remains the deciding experiment. Note the reproduction caveat: IntAct is also
literature-curated, so heavily-studied centrosome biology is over-represented there too -- the
substrate is physically-built, but not bias-free. **Test 14 maps exactly how far that caveat can
be pushed and finds a hard ceiling**: the one genuinely unbiased binary interactome (HuRI,
systematic Y2H) is structurally blind to the centrosome (10 of 14 core proteins absent, PLK1 at
degree 1), so it cannot serve as a bias-free confirmation substrate. Every dataset that covers
this biology is literature-curated; the one that isn't cannot see it. The shared-study-bias
caveat therefore cannot be removed computationally -- confirming or refuting Nsp13-PLK1 now
strictly requires a wet-lab experiment. That is the honest end of the computational road for
this hypothesis.

**Nsp4 -> TOMM70**: effectively dead. It passes the null model (not a generic hub) and negative
controls, but fails the threshold-robustness test (only appears at the loosest STRING confidence
setting), has zero external support for this specific pairing (TOMM70's real, repeatedly-confirmed
viral partner is ORF9b, already in the original 332), and now **fails to reproduce on the
independent IntAct substrate** (test 12: 0 of 8 Nsp4 preys overlap TOMM70's physical partners,
unranked). Four independent lines of evidence point the same way. This should be **withdrawn**,
not merely downgraded -- the original "compelling cross-validated finding" framing was wrong.

## The bigger revision: does the method generalize at all?

Test 10 is the most important result added after the initial validation pass, and it should
reset expectations for the whole guilt-by-association approach, not just the two headline
predictions. Test 3's internal holdout (hide 20% of Gordon et al.'s own interactions, try to
recover them) showed real signal: 11.1% recall vs 0% null. Test 10 asks a harder, more honest
question: instead of splitting one dataset into train/test folds, use Gordon et al.'s *entire*
real prey list as input and see whether the method can predict what a *separate lab, separate
cell line, separate pipeline* (Stukalov et al.) independently found. Across 18 testable baits,
mean recall dropped to 0.07% -- statistically indistinguishable from a 0.02% degree-matched
null. **The method that looked like it was working in test 3 mostly fails to generalize to
truly external biology.** The likely reason: proteins co-reported within one AP-MS study tend
to already be well-known functional partners in STRING's evidence graph (shared complexes,
shared literature, shared co-expression), so recovering a held-out member of the *same* study's
own list is an easier task than predicting a genuinely new discovery.

The one bright spot is worth naming precisely: for Orf3a, the method recovered 4 of Stukalov's
329 novel interactors -- VPS16, VPS18, VPS41, and VPS45, all HOPS complex subunits -- against
1/329 for the matched null. This lines up with real, independently published virology (ORF3a
disrupting HOPS-complex-mediated autophagosome-lysosome fusion), and it happened on the one
bait with enough independent data (329 targets) to give the test real statistical power. That's
the actual lesson: this method may work, but only demonstrably so far in the one case with
enough data to detect it. Every other bait's test was underpowered (many had fewer than 5 novel
targets to recover at all), so most of those near-zero recalls are inconclusive, not negative.

## Why the method works for some baits and not others (test 11)

Test 10 said the method mostly doesn't generalize; test 11 (`10_recall_gap_diagnosis.py`)
explains the mechanism, and it is not flattering. For each bait we pulled the STRING
sub-network *among its known preys only* (no expansion) and measured how densely those preys
are already interconnected, decomposed by evidence channel. Two things fall out:

1. **The internal-holdout recall gap is real and predictable.** Recall tracks intra-prey-set
   edge density (Spearman rho +0.65): the method recovers held-out interactions precisely for
   baits whose known interactors already form a tight STRING cluster (N, Nsp13, Nsp9), and
   recovers nothing for baits whose interactors are sparsely linked (M, Nsp12, Nsp7, Orf9b).
   Guilt-by-association can only extend a neighborhood that already exists.

2. **The clustering that drives success is literature-biased, not evidence-based.** The single
   channel most correlated with recall is *text-mining* (tscore, rho +0.78) -- co-mention in
   the literature -- while the hard-experimental-evidence channel (escore) correlates only
   +0.50, below the ~0.65 significance threshold for n=10. **Nsp8 is the decisive
   counterexample**: it has the second-highest experimental-evidence density of any bait, yet
   recovers essentially nothing in the holdout, because its interactors are not the ones the
   literature co-mentions. Baits "work" when their biology is famous and heavily co-cited, not
   when it is experimentally well-supported.

This closes the loop with test 10. A method whose successes ride on literature co-mention
density will, by construction, look good on interactions between already-famous proteins (the
internal holdout) and fail on genuinely novel, less-studied interactions (independent ground
truth) -- which is exactly the pattern observed. See `output/validation/recall_gap_diagnosis.png`.

## Honest caveats about the validation itself

- Tests 1-7 all still live inside the STRING/Enrichr/UniProt ecosystem to varying degrees --
  guilt-by-association validated by more association is a closed loop, which is exactly why
  test 8 (external, independent data) is the only one that can actually confirm rather than
  merely fail to reject. It did not confirm either headline call.
- The holdout precision numbers (test 3) are a lower bound: a "wrong" candidate could be a real,
  undiscovered interactor rather than noise -- we only have a ground truth for what the original
  screen happened to detect.
- IntAct coverage is not exhaustive; absence of an independent Nsp13-PLK1 or Nsp4-TOMM70 record
  is evidence of absence only to the extent IntAct's curation is complete for these proteins,
  which cannot be fully verified from outside.
- Test 9 underscores this from a different angle: even comparing two labs' own definitions of
  "real Nsp13 interactors" head-on, Gordon et al. (40 proteins) and Stukalov et al. (1 protein,
  not shared with Gordon's list) barely overlap at all. AP-MS interactome reproducibility across
  cell lines and labs is known to be modest in general -- a third dataset's silence on PLK1 is
  weak evidence specifically because it's also nearly silent on everything else for this bait.
  A cleaner test would need a dataset where Nsp13 produced a rich hit list, not a sparse one.
- Test 10 has the same asymmetric-power problem: 13 of 18 testable baits had fewer than 10
  novel Stukalov targets, so a 0% recall there is barely more informative than a coin flip with
  one coin. Treat the aggregate 0.07%-vs-0.02% comparison as suggestive, and the Orf3a result
  (329 targets, real signal above null) as the only individually well-powered data point.
- Test 12's positive PLK1 result is cross-*method* reproduction (STRING combined score ->
  IntAct physical interactions), not cross-*experiment* confirmation. Both graphs are still
  built from published literature, so a heavily-studied hub like PLK1 in heavily-studied
  centrosome biology is favoured in both. It rules out a STRING-specific/text-mining-only
  artifact; it does not substitute for a direct Nsp13-PLK1 binding assay.
- Test 14 carries a second methodological mismatch on top of the coverage gap: HuRI is a
  *binary* Y2H map (direct pairwise contacts), whereas the Gordon predictions concern *AP-MS
  co-complex* membership -- a protein can be in a complex without binary-binding any single
  member. Across all baits the Gordon prey sets have at most 3 internal HuRI edges, so HuRI is
  a poor substrate for this class of prediction in general, and doubly so for the centrosome,
  which it barely samples. This reinforces rather than weakens the conclusion: there is no
  existing unbiased binary resource well-matched to confirming these AP-MS-derived predictions.

## Turning the lessons into a guardrail

Two of these tests (identifier integrity, coverage) only ran because they were requested after
the fact. `pipeline/preflight.py` now wires both into a gate that runs *before* enrichment
or prediction and exits non-zero on a broken join, malformed accession, or near-total loss of
library coverage -- so any future run of this pipeline inherits the checks by default:

```bash
python3 pipeline/preflight.py && python3 pipeline/enrichment.py   # prediction gated on a clean interactome
```

## How to reproduce

```bash
source .venv/bin/activate
python3 pipeline/preflight.py                     # data-integrity gate (tests 1 & 3, automated)
python3 validation/01_identifier_integrity.py
python3 validation/02_coverage_audit.py
python3 validation/03_strip_circular.py
python3 validation/04_threshold_sweep.py
python3 validation/05_null_model.py
python3 validation/06_holdout_precision_recall.py
python3 validation/07_external_validation.py
python3 validation/08_independent_dataset_crosscheck.py
python3 validation/09_independent_ground_truth_precision_recall.py
python3 validation/10_recall_gap_diagnosis.py
python3 validation/11_alternative_substrate_intact.py
python3 validation/12_broaden_coverage.py
python3 validation/13_binary_interactome_huri.py
```

Outputs land in `data/validation/*.csv` and `output/validation/*.png`. Test 14 also writes a
reusable parsed HuRI edge list to `data/huri_edges.csv` (from NDEx UUID
`73bc2c06-5fb2-11e9-9f06-0ac135e8bacf`, owner ccsb / Vidal lab).
