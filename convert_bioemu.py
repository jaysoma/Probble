"""BioEmu sample trajectories (samples.xtc + topology.pdb) -> one multi-model PDB per protein for Probble.
Keeps up to MAXF frames, evenly spaced, so a set stays light enough to glide through in the browser.
usage: python -I convert_bioemu.py <dir with one subfolder per protein> <out dir>"""
import sys
from pathlib import Path

import mdtraj as md
import numpy as np

MAXF = 60
src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
for d in sorted(p for p in src.iterdir() if p.is_dir()):
    t = md.load(str(d / 'samples.xtc'), top=str(d / 'topology.pdb'))
    keep = np.linspace(0, t.n_frames - 1, min(MAXF, t.n_frames)).round().astype(int)
    sub = t[keep]
    name = 'SIM_' + d.name.upper()
    path = out / f'{name}.pdb'
    sub.save_pdb(str(path))           # writes MODEL/ENDMDL blocks, one per frame
    atoms = sorted({a.name for a in t.topology.atoms})
    print(f'{d.name}: {t.n_frames} samples, {t.n_residues} residues, {t.n_atoms} atoms ({", ".join(atoms[:8])}) -> {path.name}, {len(keep)} models, {path.stat().st_size/1e6:.1f} MB')
