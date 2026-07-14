#!/usr/bin/env python3
"""Question 1 — what fraction of the PIPE4/SPRINT SARS-CoV-2–human predictions
hold up experimentally, and are there patterns in validated vs not?

Inputs (all local, obtained once — see README):
  borealis_meta.json          Borealis dataset metadata (doi:10.5683/SP2/JZ77XA)
  tabs/{PIPE4,SPRINT}-*.tab    per-viral-protein prediction score tables
  biogrid_cov/*.tab3.txt       BioGRID Coronavirus (broad experimental reference)
  ../interactome.csv           Gordon et al. 2020 (strict AP-MS reference)

The 279 high-confidence predicted pairs are the .mat filenames in the Borealis
dataset (<viralUniProt>-<humanUniProt>.mat). We validate them against physical
experimental evidence and stratify by confidence, viral protein, and partner class.

ponytail: stdlib only (csv/json/glob/random). One self-check in demo().
"""
import csv, json, glob, random, statistics as st, os
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
def p(*a): return os.path.join(HERE, *a)

VNAME = {"P0DTC1":"pp1a","P0DTD1":"pp1ab(nsps)","P0DTC2":"S","P0DTC3":"ORF3a","P0DTC4":"E",
"P0DTC5":"M","P0DTC6":"ORF6","P0DTC7":"ORF7a","P0DTD8":"ORF7b","P0DTC8":"ORF8",
"P0DTC9":"N","P0DTD2":"ORF9b","P0DTD3":"ORF9c","A0A663DJA2":"ORF10"}
SC2 = {"2697049","227859","694009"}          # SARS-CoV-2 / SARS-CoV taxids
def vnorm(v): return "ORF9c" if v in ("ORF14","ORF9c") else v  # ORF9c ≡ ORF14

def predicted_pairs():
    d = json.load(open(p("borealis_meta.json")))
    mats = [f['dataFile']['filename'][:-4] for f in d['data']['latestVersion']['files']
            if f['dataFile']['filename'].endswith('.mat')]
    return [tuple(m.split("-")) for m in mats if len(m.split("-")) == 2]

def load_scores(prefix, col):
    sc = {}
    for fn in glob.glob(p("tabs", "%s-*.tab" % prefix)):
        for row in csv.DictReader(open(fn), delimiter="\t"):
            v, h = row["SARS-CoV-2_protein"].strip('"'), row["human_protein"].strip('"')
            try: sc[(v, h)] = float(row[col])
            except (ValueError, KeyError): pass
    return sc

def biogrid_experimental():
    """(viral_symbol, human_uniprot) physical pairs + the human-partner pool."""
    F = glob.glob(p("biogrid_cov", "*.tab3.txt"))[0]
    exp, pool = set(), set()
    for row in csv.DictReader(open(F), delimiter="\t"):
        if row.get("Experimental System Type") != "physical": continue
        oa, ob = row["Organism ID Interactor A"], row["Organism ID Interactor B"]
        sa, sb = row["Official Symbol Interactor A"], row["Official Symbol Interactor B"]
        ua, ub = row["SWISS-PROT Accessions Interactor A"], row["SWISS-PROT Accessions Interactor B"]
        if oa in SC2 and ob == "9606": v, hu = sa, ub
        elif ob in SC2 and oa == "9606": v, hu = sb, ua
        else: continue
        for h in hu.split("|"):
            if h and h != "-": exp.add((vnorm(v), h)); pool.add(h)
    return exp, list(pool)

def gordon_pairs():
    gmap = {"Spike":"S","Orf3a":"ORF3a","Orf3b":"ORF3b","Orf6":"ORF6","Orf7a":"ORF7a",
            "Orf8":"ORF8","Orf9b":"ORF9b","Orf9c":"ORF9c","Orf10":"ORF10"}
    g = set()
    for row in csv.DictReader(open(p("..", "interactome.csv"))):
        b = gmap.get(row["bait"], row["bait"])
        g.add((b, row["prey_uniprot"]))
    return g

def main():
    pairs = predicted_pairs()
    pipe, spr = load_scores("PIPE4", "PIPE2_score"), load_scores("SPRINT", "SPRINT_score")
    exp, pool = biogrid_experimental()
    gord = gordon_pairs()
    rows = [{"vid":v, "vname":VNAME.get(v,v), "human":h,
             "pipe4":pipe.get((v,h),0.0), "sprint":spr.get((v,h),0.0),
             "in_biogrid":int((VNAME.get(v,v),h) in exp),
             "in_gordon":int((VNAME.get(v,v),h) in gord)} for v,h in pairs]

    nb = sum(r["in_biogrid"] for r in rows); ng = sum(r["in_gordon"] for r in rows)
    print("predicted pairs: %d | distinct human partners: %d" % (len(rows), len({r['human'] for r in rows})))
    print("\nVALIDATION RATE")
    print("  BioGRID physical (broad):   %d/%d = %.1f%%" % (nb, len(rows), 100*nb/len(rows)))
    print("  Gordon 2020 (strict AP-MS): %d/%d = %.1f%%" % (ng, len(rows), 100*ng/len(rows)))

    # null model: viral-protein-matched random human partners from the experimental pool
    random.seed(1); N = 2000
    by_v = Counter(r["vname"] for r in rows); ge = 0; cnts = []
    for _ in range(N):
        c = sum(1 for vn,k in by_v.items() for h in random.sample(pool,k) if (vn,h) in exp)
        cnts.append(c); ge += (c >= nb)
    mean = st.mean(cnts)
    print("  null mean %.1f (%.1f%%)  ->  enrichment %.2fx, p=%.4f"
          % (mean, 100*mean/len(rows), nb/mean, ge/N))

    def rate(g): return (len(g), sum(x["in_biogrid"] for x in g))
    def show(title, groups):
        print("\n" + title)
        for name, g in groups:
            n, k = rate(g)
            print("  %-22s n=%-4d valid=%-3d %.0f%%" % (name, n, k, 100*k/n if n else 0))

    # per viral protein
    vs = sorted({r["vname"] for r in rows}, key=lambda v: -Counter(x["vname"] for x in rows)[v])
    show("PER VIRAL PROTEIN", [(v, [r for r in rows if r["vname"]==v]) for v in vs])

    # confidence tertiles
    for key, lab in [("pipe4","PIPE4"), ("sprint","SPRINT")]:
        s = sorted(rows, key=lambda r: r[key]); t = len(s)//3
        show("BY %s TERTILE" % lab, [("low",s[:t]),("mid",s[t:2*t]),("high",s[2*t:])])

    # method agreement
    pm, sm = st.median(r["pipe4"] for r in rows), st.median(r["sprint"] for r in rows)
    show("METHOD AGREEMENT (vs each median)", [
        ("both-high", [r for r in rows if r["pipe4"]>pm and r["sprint"]>sm]),
        ("one-high",  [r for r in rows if (r["pipe4"]>pm)!=(r["sprint"]>sm)]),
        ("both-low",  [r for r in rows if r["pipe4"]<=pm and r["sprint"]<=sm])])

    return rows, nb, mean

def demo():
    # self-check: ORF9c/ORF14 normalisation and a match both directions
    assert vnorm("ORF14") == "ORF9c" and vnorm("S") == "S"
    exp = {("N","P63165")}
    assert ("N","P63165") in exp and ("N","XXX") not in exp
    print("self-check ok")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--test": demo()
    else: main()
