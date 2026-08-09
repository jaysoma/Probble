#!/usr/bin/env python3
# Derive a CATH-style fold class for each cached AlphaFold structure from its own Ca geometry.
# Cheap secondary-structure proxy: helix => Ca(i)-Ca(i+4) is tight (~6.2A); strand => chain is locally extended.
import os, json, math, glob

DIR = os.path.dirname(os.path.abspath(__file__))
STRUCT = os.path.join(DIR, 'structures')

def cas(path):
    out = []
    with open(path) as f:
        for ln in f:
            if ln.startswith('ATOM') and ln[12:16].strip() == 'CA':
                try: out.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
                except ValueError: pass
    return out

def d(a, b): return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def classify(c):
    n = len(c)
    if n < 8: return ('Few SS', 0.0, 0.0)
    h = s = 0
    for i in range(n):
        if i+4 < n and d(c[i], c[i+4]) < 6.5: h += 1          # helix: i..i+4 close
        if i+3 < n and d(c[i], c[i+3]) > 9.5 and d(c[i], c[i+2]) > 6.0: s += 1  # strand: extended
    ph, ps = h/n, s/n
    if ph < 0.12 and ps < 0.10: cls = 'Few SS'
    elif ph >= 1.6*ps:          cls = 'Mainly alpha'
    elif ps >= 1.3*ph:          cls = 'Mainly beta'
    else:                       cls = 'Alpha + beta'
    return (cls, round(ph, 3), round(ps, 3))

names = {}
idx = os.path.join(DIR, 'index.tsv')
if os.path.exists(idx):
    with open(idx, encoding='utf-8') as f:
        next(f, None)
        for ln in f:
            p = ln.rstrip('\n').split('\t')
            if len(p) >= 3: names[p[0]] = (p[1], p[2])

rows = []
for path in sorted(glob.glob(os.path.join(STRUCT, '*.pdb'))):
    acc = os.path.splitext(os.path.basename(path))[0]
    c = cas(path)
    if not c: continue
    cls, ph, ps = classify(c)
    nm, ln = names.get(acc, (acc, str(len(c))))
    rows.append({'acc': acc, 'name': nm, 'len': len(c), 'cls': cls, 'ph': ph, 'ps': ps})

order = {'Mainly alpha': 0, 'Mainly beta': 1, 'Alpha + beta': 2, 'Few SS': 3}
rows.sort(key=lambda r: (order.get(r['cls'], 9), r['name'].lower()))

with open(os.path.join(DIR, 'catalog.json'), 'w', encoding='utf-8') as f:
    json.dump(rows, f)

from collections import Counter
dist = Counter(r['cls'] for r in rows)
print('catalog entries:', len(rows))
for k in ['Mainly alpha', 'Mainly beta', 'Alpha + beta', 'Few SS']:
    print(f'  {k:14s}: {dist.get(k,0)}')
