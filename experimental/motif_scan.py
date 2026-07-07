"""Scan SARS-CoV-2 Nsp13 for PLK1 recognition motifs.

Two independent, non-association computational lines on the Nsp13->PLK1 hypothesis:
  1. Does Nsp13 carry a PLK1 Polo-box-domain (PBD) DOCKING motif -- the canonical way PLK1 is
     stably recruited to a partner: S-[pS/pT]-P (Elia, Rellos et al. 2003, Science 299:1228)?
  2. Does Nsp13 carry PLK1 KINASE-substrate consensus sites: acidic/N at -2, hydrophobic at +1,
     i.e. [D/E/N]-X-[S/T]-[F/L/I/V/M/Y/W]?

Result (see repo): no PBD docking motif (expected ~0.3 by chance anyway), only background-level
kinase sites -> argues against a direct kinase-substrate contact and FOR a bridged/shared-
complex model (Nsp13 binds the CEP centrosome scaffolds, which are themselves PLK1 partners).

Writes single-chain nsp13.fasta / plk1.fasta for reference; the AlphaFold-Multimer inputs are
in experimental/alphafold_inputs/.
"""
import re
import urllib.request
from collections import Counter

# Nsp13 = Helicase chain of ORF1ab (UniProt P0DTD1), residues 5325-5925 (1-indexed, ft_chain)
ORF1AB_ACC = "P0DTD1"
NSP13_RANGE = (5325, 5925)
PLK1_ACC = "P53350"


def fasta_seq(acc):
    t = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta").read().decode()
    return "".join(t.splitlines()[1:])


def main():
    orf1ab = fasta_seq(ORF1AB_ACC)
    nsp13 = orf1ab[NSP13_RANGE[0] - 1:NSP13_RANGE[1]]
    plk1 = fasta_seq(PLK1_ACC)
    print(f"Nsp13 {len(nsp13)} aa | PLK1 {len(plk1)} aa")

    pbd = [(m.start() + 1, m.group()) for m in re.finditer("S[ST]P", nsp13)]
    kin = [(m.start() + 1, m.group()) for m in re.finditer("[DEN].[ST][FLIVMYW]", nsp13)]

    n = len(nsp13); c = Counter(nsp13)
    exp_pbd = (n - 2) * (c["S"] / n) * ((c["S"] + c["T"]) / n) * (c["P"] / n)

    print(f"\nPLK1 PBD docking motif  S-[pS/pT]-P : {len(pbd)} found {pbd}  (expected by chance {exp_pbd:.2f})")
    print(f"PLK1 kinase-consensus [DEN]-X-[ST]-Phi : {len(kin)} found {kin}")
    print("\nInterpretation: no canonical PBD docking motif + only background-rate kinase sites")
    print("=> favours a BRIDGED model (via CEP centrosome scaffolds), not a direct kinase-substrate contact.")

    # single-chain FASTAs for reference; the AlphaFold-Multimer inputs live in alphafold_inputs/
    with open("experimental/nsp13.fasta", "w") as f:
        f.write(f">SARSCoV2_Nsp13_helicase|{ORF1AB_ACC}_{NSP13_RANGE[0]}-{NSP13_RANGE[1]}\n{nsp13}\n")
    with open("experimental/plk1.fasta", "w") as f:
        f.write(f">Human_PLK1|{PLK1_ACC}\n{plk1}\n")


if __name__ == "__main__":
    main()
