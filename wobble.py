#!/usr/bin/env python3
"""
wobble.py -- precompute a protein's natural sway directions for the Wobble tab.

    .venv-wobble\Scripts\python wobble.py P09683 O43237        # a few accessions
    .venv-wobble\Scripts\python wobble.py --all                # every structure in structures/

For each structure this runs an Anisotropic Network Model (ANM, Bahar lab / ProDy): every
residue's alpha-carbon becomes a node, nodes within CUTOFF angstroms are joined by identical
springs, and the eigenvectors of that spring network's Hessian are the directions the
structure sways along most easily. The lowest few are the large, slow, collective motions --
domain hinges and breathing -- which is the motion a structural biologist recognises as real.

The eigen-decomposition is the expensive part (minutes for the largest chains), so it runs
here, once, and the browser only ever multiplies the stored vectors by a sine wave. The app
stays static; nothing at runtime needs Python.

Output: wobble/<ACC>.bin, little-endian --
    uint32  nAtoms                 every atom in file order, matching 3Dmol's parse
    uint32  nModes
    float32 eigval[nModes]         stiffness of each mode (softer = larger thermal amplitude)
    float32 first[3], last[3]      native x,y,z of atom 0 and atom nAtoms-1 -- the browser checks
                                   these against its own parse so a misaligned file refuses
                                   rather than jittering the wrong atoms
    float32 vec[nModes][nAtoms][3] displacement direction per atom, per mode; each mode scaled so
                                   its longest arrow is 1.0

Atoms in residues with no alpha-carbon (none in AlphaFold output) get a zero vector, so the
file always covers every atom 3Dmol will draw.
"""
import struct
import sys
import time
from pathlib import Path

import numpy as np
import prody

prody.confProDy(verbosity="none")

ROOT = Path(__file__).resolve().parent
STRUCT = ROOT / "structures"
OUT = ROOT / "wobble"
CUTOFF = 15.0       # angstrom; the ANM default and the value the method was validated at
N_MODES = 6         # modes shipped to the browser
N_CANDIDATES = 20   # modes computed before the localisation filter picks N_MODES of them
CORE_PLDDT = 70     # residues below this confidence are excluded from the spring network


def pdb_path(acc: str) -> Path:
    # the reference protein lives at the project root under a different name
    if acc == "O05011":
        return ROOT / "O05011_af.pdb"
    return STRUCT / f"{acc}.pdb"


def compute(acc: str) -> Path:
    ag = prody.parsePDB(str(pdb_path(acc)))
    n_all = ag.numAtoms()
    betas = ag.getBetas()                         # pLDDT rides the B-factor column

    # The network is built on the CONFIDENT core only. AlphaFold's low-confidence tails are
    # extended chains with two neighbours each; in the spring network they bend for nearly
    # nothing, so they would fill every soft mode with a tail flapping while the body sat
    # still. Excluding them leaves the modes describing the fold itself.
    core = ag.select(f"calpha and beta >= {CORE_PLDDT}")
    if core is None or core.numAtoms() < 3:
        raise RuntimeError(f"{acc}: too few confident alpha-carbons")

    anm = prody.ANM(acc)
    anm.buildHessian(core, cutoff=CUTOFF)
    anm.calcModes(n_modes=N_CANDIDATES)

    # Localisation filter: even the core can carry a loose loop. Score each candidate by how
    # much of its motion sits in its most-moving 5% of nodes and keep the most spread-out
    # N_MODES, softest first among those.
    V = anm.getEigvecs().T.reshape(N_CANDIDATES, -1, 3)
    E = (V ** 2).sum(2); E = E / E.sum(1, keepdims=True)
    top5 = np.array([np.sort(e)[::-1][: max(1, e.size // 20)].sum() for e in E])
    keep = np.argsort(top5)[:N_MODES]
    keep = keep[np.argsort(anm.getEigvals()[keep])]
    sub = anm[keep.tolist()]

    # Carry each confident residue's node vector out to all its atoms, rigidly.
    ext, ext_atoms = prody.extendModel(sub, core, ag.select(f"beta >= {CORE_PLDDT}"))
    full = np.zeros((N_MODES, n_all, 3), dtype=np.float32)
    full[:, ext_atoms.getIndices(), :] = ext.getEigvecs().T.reshape(N_MODES, -1, 3)

    # Tail atoms ride the nearest confident alpha-carbon's vector: they follow the body's sway
    # instead of dictating a fake one. The browser adds their own uncertainty cloud on top.
    low = np.where(betas < CORE_PLDDT)[0]
    if low.size:
        xyz, cxyz = ag.getCoords(), core.getCoords()
        nearest = ((xyz[low, None, :] - cxyz[None, :, :]) ** 2).sum(2).argmin(1)
        node_vec = sub.getEigvecs().T.reshape(N_MODES, -1, 3)
        full[:, low, :] = node_vec[:, nearest, :]

    # Per mode, normalise so the longest arrow is 1.0; the browser sets amplitude in angstrom.
    for m in range(N_MODES):
        full[m] /= max(1e-9, np.linalg.norm(full[m], axis=1).max())

    coords = ag.getCoords().astype(np.float32)
    eigvals = anm.getEigvals()[keep].astype(np.float32)

    OUT.mkdir(exist_ok=True)
    dst = OUT / f"{acc}.bin"
    with open(dst, "wb") as f:
        f.write(struct.pack("<II", n_all, N_MODES))
        f.write(eigvals.tobytes())
        f.write(coords[0].tobytes())
        f.write(coords[-1].tobytes())
        f.write(np.ascontiguousarray(full).tobytes())
    return dst


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    accs = ([p.stem for p in sorted(STRUCT.glob("*.pdb"))] + ["O05011"]) if argv == ["--all"] else argv
    failed = 0
    for i, acc in enumerate(accs, 1):
        t0 = time.perf_counter()
        try:
            dst = compute(acc.upper())
            print(f"[{i}/{len(accs)}] {acc:8s} {dst.stat().st_size/1e6:6.2f} MB  {time.perf_counter()-t0:6.1f}s", flush=True)
        except Exception as e:  # keep the batch going; report at the end
            failed += 1
            print(f"[{i}/{len(accs)}] {acc:8s} FAILED: {e}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
