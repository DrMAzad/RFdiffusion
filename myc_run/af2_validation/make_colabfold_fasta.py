#!/usr/bin/env python
"""Build ColabFold complex FASTAs (target:binder) from ProteinMPNN output.

Usage: python make_colabfold_fasta.py <mpnn_seqs_dir> <out.fasta> [--all]
  Default: one entry per design using its best (lowest global_score) sequence.
  --all:   every ProteinMPNN sample (name_sN), for deeper screening.
"""
import sys, os, re, glob

# c-Myc leucine-zipper target used for design (1NKP chain A, residues 955-984;
# C-terminal residue is the reactive Cys984). Override with env MYC_TARGET_SEQ.
TARGET = os.environ.get("MYC_TARGET_SEQ", "AEEQKLISEEDLLRKRREQLKHKLEQLGGC")

def records(fa):
    recs, cur = [], None
    for ln in open(fa):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if cur: recs.append(cur)
            cur = {"h": ln[1:], "s": ""}
        elif cur and ln.strip():
            cur["s"] += ln.strip()
    if cur: recs.append(cur)
    out = []
    for r in recs:
        if "sample=" not in r["h"]:
            continue
        m = re.search(r"global_score=([0-9.eE+-]+)", r["h"])
        sm = re.search(r"sample=(\d+)", r["h"])
        out.append({"gs": float(m.group(1)) if m else 9e9,
                    "sample": int(sm.group(1)) if sm else 0,
                    "binder": min(r["s"].split("/"), key=len)})
    return out

def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    seqs_dir, out = sys.argv[1], sys.argv[2]
    want_all = "--all" in sys.argv[3:]
    n = 0
    with open(out, "w") as fh:
        for fa in sorted(glob.glob(os.path.join(seqs_dir, "*.fa"))):
            name = os.path.basename(fa)[:-3]
            recs = records(fa)
            if not recs:
                continue
            chosen = recs if want_all else [min(recs, key=lambda r: r["gs"])]
            for r in chosen:
                tag = name if not want_all else "%s_s%d" % (name, r["sample"])
                fh.write(">%s\n%s:%s\n" % (tag, TARGET, r["binder"]))
                n += 1
    print("wrote %d complex entries -> %s  (target len %d)" % (n, out, len(TARGET)))

if __name__ == "__main__":
    main()
