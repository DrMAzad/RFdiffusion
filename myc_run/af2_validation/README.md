# AlphaFold validation of Myc binder designs (GPU step)

The CPU pipeline (RFdiffusion → ProteinMPNN → MM interaction energy) *proposes and
ranks* candidates. The decisive in-silico filter before synthesis is **AlphaFold**:
predict each binder:Myc complex and keep only designs AF2 is confident about. This
needs a GPU — run it on a rented cloud GPU, an HPC node, or Google Colab.

## The metric that matters: `pae_interaction`

For RFdiffusion binders (Bennett et al., *Nat. Commun.* 2023), the single best
computational predictor of real binding is the **interface predicted aligned error** —
the mean AF2 PAE between binder residues and target residues (both directions). Canonical
filters:

| Metric | Keep if | Meaning |
|---|---|---|
| **pae_interaction** | **< 10** (ideally < 5) | AF2 is confident the two chains sit together as modeled |
| binder pLDDT | > 80 | the peptide is well-folded, not floppy |
| ipTM (multimer) | > 0.5 | overall complex confidence |
| binder RMSD to design | < ~2 Å | AF2 agrees with the RFdiffusion backbone |

For this **covalent** project, also carry forward the CPU metric: the design must still
place a residue within warhead reach (~3–4 Å) of **Cys984-SG** (see `../binding_assessment.csv`).

## Two routes

### Route A — ColabFold (easiest, no license)
`localcolabfold` on any GPU, or the free ColabFold notebook. Predicts the
binder:target complex from sequence.

```bash
# 1. build target:binder complex FASTAs from the ProteinMPNN sequences
python make_colabfold_fasta.py ../mpnn/seqs complexes.fasta
# 2. predict (GPU). Install: https://github.com/YoshitakaMo/localcolabfold
colabfold_batch --num-recycle 3 --model-type alphafold2_multimer_v3 complexes.fasta af2_out
# 3. rank by pae_interaction + binder pLDDT
python parse_af2.py af2_out af2_metrics.csv
```

### Route B — AF2 initial-guess (gold standard for binders)
The Baker-lab `dl_binder_design` pipeline seeds AF2 with the RFdiffusion backbone, which
is markedly more accurate for binders than vanilla AF2. Heavier setup (AF2 params + the
repo). See `run_af2_initial_guess.sh` for the exact steps; it consumes the RFdiffusion
PDBs + ProteinMPNN sequences directly and writes a scorefile whose `pae_interaction`
column you filter on the same way.

## Reading the result
`parse_af2.py` writes a CSV ranked by `pae_interaction` (lower = better). Synthesize the
designs that pass **all** filters AND keep a short Cys984-SG distance. Expect only a
handful of the scan to survive — that is normal and the point of the filter.
