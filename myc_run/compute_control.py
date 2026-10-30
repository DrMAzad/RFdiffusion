#!/usr/bin/env python
"""Positive controls: score the paper's MAX-derived covalent binder (compound 24/26,
the MYC-binding entity of compounds 27 & 28) and native MAX, in the native MYC-MAX
binding pose from 1NKP, with the SAME MM pipeline used for the AI peptides.

Binder is built by threading the control sequence onto the MAX chain-B backbone
(C-terminal 24 residues, which pack against MYC Cys438/Cys984). The non-natural
Dap(N-acrylamide) warhead X is modeled as Ala (its Cbeta marks the warhead stem;
the acrylamide extends ~3 A further toward the thiol). Target = MYC chain A 955-984,
with native sidechains kept."""
import os, tempfile, numpy as np
from pdbfixer import PDBFixer
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile
from openmm import LangevinIntegrator, Context, Platform, LocalEnergyMinimizer
from openmm.unit import kilocalorie_per_mole, nanometer, kelvin, picosecond, femtosecond
import freesasa

NKP = "/home/user/myc_struct/1NKP.pdb"
MYC_LO, MYC_HI, CYS = 955, 984, 984            # MYC chain A target window + reactive Cys
MAX_CTERM_N = 24                                # MAX chain B C-terminal residues used as scaffold
AA1TO3 = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY',
          'H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER',
          'T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
# compound 24 = clean 24-mer MAX-derived covalent binder (X->A); 26/27/28 share this binding seq
CTRLS = {
    "cmpd24_26 (paper)": "EIDDLKREIALLEQEIRALEKAAA",   # X(Dap-acrylamide) modeled as A
    "native_MAX":        None,                          # keep native MAX sequence/atoms
}
ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
plat = Platform.getPlatformByName("CPU")
freesasa.setVerbosity(freesasa.silent)

def read_atoms():
    myc, maxb = [], []
    for ln in open(NKP):
        if not ln.startswith("ATOM"): continue
        ch, ri, alt = ln[21], int(ln[22:26]), ln[16]
        if alt not in (" ", "A"): continue
        if ch == "A" and MYC_LO <= ri <= MYC_HI: myc.append(ln)
        elif ch == "B": maxb.append(ln)
    return myc, maxb

def build_pdb(myc_lines, maxb_lines, seq, out):
    # MAX B: take C-terminal N residues, backbone only, renumber 1..N, chain B
    resis = sorted({int(l[22:26]) for l in maxb_lines})[-MAX_CTERM_N:]
    keep_bb = {"N", "CA", "C", "O"}
    with open(out, "w") as fh:
        for l in myc_lines:                      # target: native, full atoms, chain A
            fh.write(l)
        fh.write("TER\n")
        serial = 1
        for newi, ri in enumerate(resis, start=1):
            rlines = [l for l in maxb_lines if int(l[22:26]) == ri and l[12:16].strip() in keep_bb]
            if seq is None:
                resn = rlines[0][17:20] if rlines else "GLY"
            else:
                resn = AA1TO3[seq[newi-1]]
            for l in rlines:
                fh.write("ATOM  %5d %s %s B%4d%s" % (serial, l[12:16], resn, newi, l[26:]))
                serial += 1
        fh.write("TER\nEND\n")

def energy_ctx(topology, positions):
    system = ff.createSystem(topology, nonbondedMethod=NoCutoff, constraints=None)
    integ = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtosecond)
    return Context(system, integ, plat)

def sasa(path):
    return freesasa.calc(freesasa.Structure(path)).totalArea()

def score(label, seq, myc_lines, maxb_lines):
    tmp = tempfile.mkdtemp()
    raw = os.path.join(tmp, "raw.pdb")
    build_pdb(myc_lines, maxb_lines, seq, raw)
    fixer = PDBFixer(filename=raw)
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
    mod = Modeller(fixer.topology, fixer.positions)
    ctx = energy_ctx(mod.topology, mod.positions); ctx.setPositions(mod.positions)
    LocalEnergyMinimizer.minimize(ctx, maxIterations=400)
    st = ctx.getState(getPositions=True, getEnergy=True)
    Ec = st.getPotentialEnergy().value_in_unit(kilocalorie_per_mole); pos = st.getPositions()
    def chainE(keep):
        m = Modeller(fixer.topology, pos)
        m.delete([a for a in m.topology.atoms() if a.residue.chain.id != keep])
        c = energy_ctx(m.topology, m.positions); c.setPositions(m.positions)
        return c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie_per_mole), m
    Eb, mb = chainE("B"); Et, mt = chainE("A")
    dE = Ec - Eb - Et
    pc, pb, pt = (os.path.join(tmp, x) for x in ("c.pdb","b.pdb","t.pdb"))
    for p,(top,ps) in zip((pc,pb,pt),((fixer.topology,pos),(mb.topology,mb.positions),(mt.topology,mt.positions))):
        with open(p,"w") as fh: PDBFile.writeFile(top, ps, fh)
    bsa = max(0.0,(sasa(pb)+sasa(pt)-sasa(pc))/2.0)
    bp,tp,sg=[],[],None
    for a in fixer.topology.atoms():
        if a.element is None or a.element.symbol=="H": continue
        p=np.array(pos[a.index].value_in_unit(nanometer))*10
        if a.residue.chain.id=="B": bp.append(p)
        else:
            tp.append(p)
            if a.residue.id==str(CYS) and a.name=="SG": sg=p
    bp,tp=np.array(bp),np.array(tp)
    dmat=np.sqrt(((bp[:,None,:]-tp[None,:,:])**2).sum(-1))
    contacts=int((dmat<4.5).sum())
    sgd=float(np.sqrt(((bp-sg)**2).sum(-1)).min()) if sg is not None else float("nan")
    L = len(seq) if seq else MAX_CTERM_N
    print("%-20s %3d %-26s %8.1f %8.0f %9d %10.1f" % (label, L, seq or "(native MAX)", dE, bsa, contacts, sgd))
    return dict(name=label,L=L,seq=seq or "native_MAX",dE=dE,bsa=bsa,contacts=contacts,sg=sgd)

def main():
    myc, maxb = read_atoms()
    print("%-20s %3s %-26s %8s %8s %9s %10s" % ("control","L","sequence","dE_kcal","BSA_A2","contacts","Cys984_A"))
    print("-"*95)
    rows=[score(k,v,myc,maxb) for k,v in CTRLS.items()]
    with open("/home/user/myc_run/control_assessment.csv","w") as fh:
        fh.write("control,length,sequence,dE_interaction_kcal_mol,buried_SASA_A2,interface_contacts,min_dist_to_Cys984_SG_A\n")
        for r in rows:
            fh.write("%s,%d,%s,%.2f,%.1f,%d,%.2f\n"%(r["name"],r["L"],r["seq"],r["dE"],r["bsa"],r["contacts"],r["sg"]))
    print("\nCSV -> /home/user/myc_run/control_assessment.csv")

if __name__=="__main__":
    main()
