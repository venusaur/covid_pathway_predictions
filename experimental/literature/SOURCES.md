# Literature dive — source log

Every source examined for the Nsp13 → PLK1 / centrosome open questions, with a relevance
verdict. Verdicts map to the hypothesis "Nsp13 engages the PLK1-regulated centrosome":
CONFIRMS / SUPPORTS / NEUTRAL / REFUTES, or SOURCE (primary dataset) / METHOD (background).
Synthesis and per-question leads: `LITERATURE_REVIEW.md`.

Note on provenance: the automated deep-research harness (run wf_7ec7e38d-16e) hit a session
rate limit and completed only 8 of 40 agents, so this table is a focused manual sweep, not the
exhaustive adversarially-verified pass. Re-run the harness after the reset for full verification.

| # | Source (short) | URL / ID | Q | Verdict | Note |
|---|---|---|---|---|---|
| 1 | Gordon et al. 2020, *Nature* — SARS-CoV-2 interactome | PMC7239059 | 1,4,5 | SOURCE | The AP-MS the whole project builds on. Nsp13 binds centrosome (CEP250, CEP68, AKAP9) + kinases TBK1, TBKBP1, PRKACA. No PLK1 |
| 2 | Stukalov et al. 2021, *Nature* — multilevel proteomics | doi:10.1038/s41586-021-03493-4 | 1 | NEUTRAL | Independent AP-MS; Nsp13 recovered only RNF8 (sparse, low power). No PLK1 |
| 3 | Chen et al. 2021, *Biochem. Genet.* — gut DEG/hub genes | PMC8596852 | 3 | NEUTRAL (weak) | PLK1 an upregulated "hub gene" — the study-bias artifact, not evidence. Expression != interaction != requirement |
| 4 | NSP13 hijacks EWSR1 to boost RNA unwinding | PMC8868009 | 1,5 | NEUTRAL | Primary, experimental. A real functional Nsp13 host partner (EWSR1), unrelated to PLK1/centrosome |
| 5 | SARS-CoV NSP13 + DNA pol delta (p125/p50) -> S-phase arrest | PMC9937006 | 2,5 | SUPPORTS (adjacent) | SARS-CoV ortholog drives S-phase arrest via Pol delta, synergizing with NSP12. Cell-cycle engagement, but the DNA-replication/S-phase axis, not PLK1/G2/M |
| 6 | P53-independent G1 arrest increases SARS-CoV-2 replication | PMC10972258 | 2 | SUPPORTS (Q2 only) | Virus perturbs cell cycle; G1 arrest aids replication. No PLK1/Nsp13 mechanism |
| 7 | Host cell-cycle checkpoint as antiviral target, *STTT* 2022 | doi:10.1038/s41392-022-01296-1 | 2,3 | SUPPORTS (Q2) | SARS-CoV-2 induces S and G2/M arrest; cell cycle as antiviral target. Doesn't name PLK1/Nsp13 |
| 8 | N protein delays cell cycle in S-phase | S143842212500027X | 2 | SUPPORTS (Q2) | Cell-cycle disruption driven by N, not Nsp13 |
| 9 | Proximal proteome (BioID) of 17 SARS-CoV-2 proteins, May et al. | PMC8513853 / PMC7924263 | 1,4 | NEUTRAL (N/A) | Checked directly: **Nsp13 is NOT one of the 17 baits**, so it says nothing about Nsp13. My earlier "high-value lead" flag was wrong |
| 14 | NSP13 binds STAT1, suppresses IFN (blocks JAK1->STAT1), co-IP + IF | PMC8574307 / PMC8939493 | 1,5 | NEUTRAL | Independent, experimentally-validated Nsp13 host partner (STAT1). Innate immunity, not PLK1/centrosome |
| 15 | PKA-CREB1 axis regulates coronavirus via Nsp13 association | PMC11019953 | 4,5 | SUPPORTS (adjacent) | Experimentally validates Nsp13-PKA(PRKACA) — a Gordon-predicted Nsp13 kinase interactor — and its functional role. PKA is centrosome-anchored (via AKAP9), so any centrosome link may run through validated PKA, not PLK1 |
| 16 | NSP13 interacts with HSP70/HSP90 and ODC (polyamine) in replication | PMC12840268 | 1,5 | NEUTRAL | More validated Nsp13 partners (chaperones, polyamine). None are PLK1/centrosome |
| 10 | Subcellular localization of SARS-CoV-2 proteins, *STTT* 2020 | doi:10.1038/s41392-020-00372-8 | 4 | NEUTRAL | Nsp13 nuclear + cytoplasmic, colocalizes with SC35 (splicing speckles). No centrosome IF reported |
| 11 | BI-2536 / BI-6727 reactivate latent HIV (dual PLK/bromodomain) | doi:10.1038/s41598-018-21942-5 | 3 | REFUTES (drug angle) | PLK1 inhibitors' antiviral use is HIV latency, not SARS-CoV-2. No SARS-CoV-2 antiviral data for BI2536/volasertib |
| 12 | Landscape of human cancer proteins targeted by SARS-CoV-2, *Cancer Discov.* 2020 | aacrjournals cancerdiscovery 10/7/916 | 3 | NEUTRAL | Reanalysis of Gordon data through a cancer/cell-cycle lens; derivative, no new PLK1 interaction data |
| 13 | Role of the early secretory pathway in SARS-CoV-2, *JCB* 2020 | rupress e202006005 | 4 | NEUTRAL | Nsp13's Golgi/centrosome links read as trafficking hijack, not mitotic |
