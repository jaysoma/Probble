"""All of a BioEmu set's shapes (10,000), packed for the Ensemble tab's build mode (Jay, 2026-11-25: "feed all 10000").
The build only needs which 1 A cell each atom lands in, so each atom of each shape is stored as three cell numbers,
one byte each, after every shape is superposed (C-alphas, then again on the steadier half). Gzipped; the browser
unpacks it with DecompressionStream.

Layout (little-endian): magic | u32 shapes | u32 atoms | f32 origin x,y,z | f32 cell | i16 resi[atoms]
                        | cell xyz[shapes][atoms][3], u8 when magic is 'PBB1', u16 when 'PBB2' (a set wider than 256 cells:
                        some BioEmu shapes are nearly fully unfolded, and adenylate kinase spans 831 A)
usage: python -I build_bins.py <folder with samples.xtc + topology.pdb> <out .build.bin.gz>"""
import gzip
import struct
import sys
from pathlib import Path

import mdtraj as md
import numpy as np

CELL = 1.0
src, out = Path(sys.argv[1]), Path(sys.argv[2])
t = md.load(str(src / 'samples.xtc'), top=str(src / 'topology.pdb'))
ca = t.topology.select('name CA')
t.superpose(t, frame=0, atom_indices=ca)
xyz = t.xyz * 10.0                                                   # nm -> angstrom
rmsf = np.sqrt(((xyz[:, ca] - xyz[:, ca].mean(0)) ** 2).sum(2).mean(0))
steady = ca[np.argsort(rmsf)[: max(3, len(ca) // 2)]]
t.superpose(t, frame=0, atom_indices=steady)
xyz = t.xyz * 10.0
lo = xyz.reshape(-1, 3).min(0)
cells = np.floor((xyz - lo) / CELL).astype(np.int64)
span = cells.max()
magic, cells = (b'PBB1', cells.astype(np.uint8)) if span <= 255 else (b'PBB2', cells.astype(np.uint16))
resi = np.array([a.residue.resSeq for a in t.topology.atoms], dtype=np.int16)
head = magic + struct.pack('<II', t.n_frames, t.n_atoms) + struct.pack('<4f', *lo, CELL) + resi.tobytes()
with gzip.open(out, 'wb', compresslevel=6) as f:
    f.write(head)
    f.write(cells.tobytes())
print(f'{src.name}: {t.n_frames} shapes x {t.n_atoms} atoms, span {span + 1} cells ({magic.decode()}) -> {out.name} {out.stat().st_size / 1e6:.1f} MB')
