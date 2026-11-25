"""Download one lab-measured ensemble (a multi-model NMR entry) per library protein, for Probble's private
Ensemble tab. For each UniProt accession, candidates are tried longest-first; an entry counts only if it holds
2 or more models. Hydrogens are stripped (3Dmol drops them anyway). Writes FoldView/ensembles/<PDB>.pdb and
FoldView/ensembles.json {accession: {pdb, models, residues}}."""
import json
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(r'C:\Projects\GithubRoot\FoldView')
OUT = ROOT / 'ensembles'
OUT.mkdir(exist_ok=True)
cands = json.load(open(HERE / 'library_nmr.json', encoding='utf-8'))     # {acc: [[pdb, sample_len], ...]}


def fetch(pid):
    req = urllib.request.Request(f'https://files.rcsb.org/download/{pid}.pdb', headers={'User-Agent': 'research-client'})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode('utf-8', 'replace')


def strip(text):
    keep, models, res = [], 0, set()
    for line in text.splitlines():
        rec = line[:6]
        if rec.startswith('MODEL'):
            models += 1
        if rec.startswith(('ATOM', 'HETATM')):
            if line[76:78].strip() == 'H' or line[12:16].strip().startswith('H') and line[76:78].strip() in ('', 'H'):
                continue
            if rec.startswith('ATOM') and models <= 1 and line[12:16].strip() == 'CA':
                res.add((line[21], line[22:27]))
        if rec.startswith(('ATOM', 'HETATM', 'MODEL', 'ENDMDL', 'TER', 'END')):
            keep.append(line)
    return '\n'.join(keep) + '\n', models, len(res)


index, tried = {}, 0
for acc, lst in sorted(cands.items()):
    seen = []
    for pid, ln in sorted(lst, key=lambda x: -x[1]):
        if pid in seen:
            continue
        seen.append(pid)
        if len(seen) > 4:
            break
        try:
            text = fetch(pid); tried += 1
        except Exception as e:
            print('  ', acc, pid, 'failed:', e); continue
        body, models, nres = strip(text)
        time.sleep(0.3)
        if models >= 2:
            (OUT / f'{pid}.pdb').write_text(body, encoding='utf-8')
            index[acc] = {'pdb': pid, 'models': models, 'residues': nres}
            print(f'{acc} -> {pid}: {models} models, {nres} residues', flush=True)
            break
    else:
        print(f'{acc}: no multi-model entry among {seen}', flush=True)
(ROOT / 'ensembles.json').write_text(json.dumps(index, indent=1), encoding='utf-8')
size = sum(f.stat().st_size for f in OUT.glob('*.pdb')) / 1e6
print(f'\n{len(index)} of {len(cands)} proteins have a usable ensemble; {tried} files fetched; {size:.0f} MB on disk')
