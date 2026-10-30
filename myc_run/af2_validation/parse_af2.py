#!/usr/bin/env python
"""Rank ColabFold predictions by pae_interaction + binder pLDDT.

Usage: python parse_af2.py <colabfold_out_dir> <out.csv>
Assumes complexes were built target-first (target:binder), as make_colabfold_fasta.py does.
pae_interaction = mean PAE across the binder<->target block of the top-ranked model.
Pass/keep thresholds follow Bennett et al. 2023: pae_interaction<10, binder pLDDT>80.
"""
import sys, os, re, json, glob

PAE_CUT, PLDDT_CUT = 10.0, 80.0

def load_scores(path):
    d = json.load(open(path))
    return d.get("pae") or d.get("predicted_aligned_error"), d.get("plddt"), d.get("ptm"), d.get("iptm")

def target_len_from_fasta(outdir):
    # ColabFold writes <name>.a3m / config; fall back to env or None
    return None

def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    outdir, out = sys.argv[1], sys.argv[2]
    tlen_env = os.environ.get("MYC_TARGET_LEN")
    rows = []
    # one top-ranked score json per design
    for sj in sorted(glob.glob(os.path.join(outdir, "*_scores_rank_001_*.json")) +
                     glob.glob(os.path.join(outdir, "*rank_001*.json"))):
        name = re.split(r"_scores_rank|_unrelaxed_rank|_rank_", os.path.basename(sj))[0]
        pae, plddt, ptm, iptm = load_scores(sj)
        if pae is None or plddt is None:
            continue
        n = len(plddt)
        # infer target length: prefer env; else the known 30-mer Myc target
        tlen = int(tlen_env) if tlen_env else 30
        if tlen >= n:
            continue
        bl = list(range(tlen, n))        # binder residue indices
        tl = list(range(0, tlen))        # target residue indices
        inter = []
        for i in bl:
            for j in tl:
                inter.append(pae[i][j]); inter.append(pae[j][i])
        pae_i = sum(inter) / len(inter)
        binder_plddt = sum(plddt[i] for i in bl) / len(bl)
        rows.append(dict(name=name, pae_interaction=pae_i, binder_plddt=binder_plddt,
                         ptm=ptm, iptm=iptm, binder_len=len(bl)))
    rows.sort(key=lambda r: r["pae_interaction"])
    with open(out, "w") as fh:
        fh.write("design,pae_interaction,binder_plddt,ptm,iptm,binder_len,PASS\n")
        for r in rows:
            passed = (r["pae_interaction"] < PAE_CUT and r["binder_plddt"] > PLDDT_CUT)
            fh.write("%s,%.2f,%.1f,%s,%s,%d,%s\n" % (
                r["name"], r["pae_interaction"], r["binder_plddt"],
                ("%.3f" % r["ptm"]) if r["ptm"] is not None else "",
                ("%.3f" % r["iptm"]) if r["iptm"] is not None else "",
                r["binder_len"], "yes" if passed else "no"))
    npass = sum(1 for r in rows if r["pae_interaction"] < PAE_CUT and r["binder_plddt"] > PLDDT_CUT)
    print("parsed %d designs, %d pass (pae_interaction<%.0f & binder pLDDT>%.0f) -> %s" %
          (len(rows), npass, PAE_CUT, PLDDT_CUT, out))

if __name__ == "__main__":
    main()
