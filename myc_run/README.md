# Short covalent peptide binder to c-Myc (Cys984)

Design a **short, linear (non-cyclic) peptide** that docks onto the c-Myc
leucine zipper and places a warhead next to the reactive cysteine — a downsized
version of the MAX-derived electrophilic helical-zipper peptides in Araújo et
al., *Angew. Chem. Int. Ed.* 2026.

## Target

- Structure: **1NKP** (c-Myc : Max : DNA). c-Myc = chains A/D, Max = B/E.
- The c-Myc chain has **exactly one cysteine, Cys984** (construct numbering
  897–984) — the reactive zipper cysteine and the covalent anchor.
- Interface knobs the peptide should engage (Myc residues that Max contacts,
  nearest the cysteine): **Leu978, Leu981, Cys984** (and Lys977 nearby).

Prepared target files (Myc chain A only, Max + DNA removed):
- `myc_target_zipperC.pdb` — C-terminal zipper, residues 955–984 (fast, focused).
- `myc_target_full.pdb` — full bHLHLZ helix, residues 897–984 (more context).

## Run

`run_myc_binder.sh` scans binder lengths 6→14 against the trimmed target with the
cysteine as hotspot, with `noise_scale=0` (the recommended setting for PPI/binder
design — cleaner backbones). From the root of an RFdiffusion checkout with a GPU,
weights in `./models`, and the `SE3nv` env active:

```bash
bash myc_run/run_myc_binder.sh
```

"As short as possible": RFdiffusion rarely forms a stable interface below ~8–10
residues. Take the **shortest** length that still passes downstream validation
rather than assuming 6 will work.

## After RFdiffusion (sequence + validation + warhead)

1. **ProteinMPNN** on each backbone to design the peptide sequence (fix nothing on
   the target).
2. **AF2 / AlphaFold-Multimer** (or ESMFold) to check the peptide refolds onto Myc —
   keep designs with high pLDDT and low interface pAE.
3. **Warhead**: RFdiffusion models no covalent bond. Pick the binder residue whose
   Cβ points at Cys984-SG, mutate it to Cys (or an electrophilic unnatural residue),
   and install the electrophile (e.g. chloroacetamide / acrylamide) so it reaches
   the Myc thiol. Minimize and check the SG–warhead distance (~3–4 Å).

## Watching it generate

Every run writes per-step trajectories to `outputs/.../traj/`:
- `*_pX0_traj.pdb` — the model's running prediction of the final structure.
- `*_Xt-1_traj.pdb` — the noisy sample at each step.

Load either multi-model file into the **Diffusion Trajectory Viewer** to watch the
peptide condense onto Myc. The generated binder is chain A (highlighted); Myc is the
muted target. Use the "reverse direction" toggle if denoising plays backwards.

## Where to run (no local GPU)

- **Google Colab** — the ColabDesign / RFdiffusion notebook gives a free (time-limited)
  GPU; good for a scan this small. Upload `myc_target_zipperC.pdb` and use the same
  contig/hotspots.
- **Rented cloud GPU** (RunPod / Lambda / Vast.ai, ~a few $/hr) — build the image from
  this repo's `docker/Dockerfile`, download weights, run the script. Best for production.
- **Institutional HPC** — Apptainer/Singularity from the same Dockerfile on a GPU node.
