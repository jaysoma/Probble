"""BioEmu sample trajectories (samples.xtc + topology.pdb) -> one multi-model PDB per set for Probble.
Finds every folder under <src> holding a samples.xtc; keeps up to MAXF shapes, evenly spaced, so a set stays light
enough to glide through in the browser. Names the output SIM_<FOLDER>.pdb (upper case, dots to underscores) and
prints one JSON line per set: {"pdb", "folder", "group", "models", "len"}.
usage: python -I convert_bioemu.py <src dir> <out dir>"""
import json
import sys
from pathlib import Path

import mdtraj as md
import numpy as np

MAXF = 60
src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
for xtc in sorted(src.rglob('samples.xtc')):
    d = xtc.parent
    t = md.load(str(xtc), top=str(d / 'topology.pdb'))
    keep = np.linspace(0, t.n_frames - 1, min(MAXF, t.n_frames)).round().astype(int)
    name = 'SIM_' + d.name.upper().replace('.', '_').replace('-', '_')
    t[keep].save_pdb(str(out / f'{name}.pdb'))
    print(json.dumps({'pdb': name, 'folder': d.name, 'group': d.parent.name, 'models': int(len(keep)), 'len': int(t.n_residues)}), flush=True)
