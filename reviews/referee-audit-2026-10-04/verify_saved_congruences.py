#!/usr/bin/env python3
"""Verify reviewer PSD certificates with no author Python imports and no floats.
Usage: python verify_saved_congruences.py /path/to/werner-problem6-proof
"""
import sys,json,math,hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
root=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/werner_review/werner-problem6-proof')
local=Path(__file__).resolve().parent
c=root/'certificate';expected=json.loads((local/'hash_check.json').read_text())['certificate_manifest_sha256']
assert hashlib.sha256((c/'manifest.json').read_bytes()).hexdigest()==expected
manifest=json.loads((c/'manifest.json').read_text())['files'];data=json.loads((c/'certificate.json').read_text())
proofs=json.loads((local/'psd_congruence_certificates.json').read_text())['proofs']
assert len(proofs)==len(data['retained_multipliers'])==716
for record,proof in zip(data['retained_multipliers'],proofs):
 path=c/record['path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==manifest[record['path']]==record['sha256']
 block=json.loads(path.read_text());assert block['id']==proof['id']==record['id']
 values=[[Fraction(x) for x in row] for row in block['Z']];d=math.lcm(*(v.denominator for row in values for v in row))
 z=np.array([[v.numerator*(d//v.denominator) for v in row] for row in values],dtype=object)
 r=np.array(proof['R'],dtype=object);n=len(z)
 assert z.shape==(n,n)==r.shape and np.array_equal(z,z.T)
 assert all(r[i,i]!=0 for i in range(n)) and all(r[i,j]==0 for i in range(n) for j in range(i+1,n))
 h=r@z@r.T;margins=[int(h[i,i])-sum(abs(int(h[i,j])) for j in range(n) if j!=i) for i in range(n)]
 assert min(margins)>0 and str(min(margins))==proof['minimum_integer_dominance_margin']
print('PASS: all 716 saved congruences establish exact positive definiteness; no floating arithmetic used.')
