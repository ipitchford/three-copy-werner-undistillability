"""Recheck the supplied 716 positive-definiteness witnesses against final data.

No imports from the certificate implementation and no floating arithmetic.
This producer check does not authenticate the external reviewer's identity.
"""
from fractions import Fraction
from pathlib import Path
import hashlib,json,math,sys
import numpy as np
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
root=Path(__file__).resolve().parent
def require(test,message):
    if not test:raise RuntimeError(message)
catalogue=json.loads((root/'certificate/certificate.json').read_text())['retained_multipliers']
audit=root/'reviews/referee-audit-2026-10-04'
proofs=json.loads((audit/'psd_congruence_certificates.json').read_text())['proofs']
require(len(catalogue)==len(proofs)==716,'Complete witness count')
for rec,proof in zip(catalogue,proofs):
    path=root/'certificate'/rec['path']
    require(hashlib.sha256(path.read_bytes()).hexdigest()==rec['sha256'],'Multiplier input binding')
    z=json.loads(path.read_text());require(z['id']==proof['id'],'Witness identity')
    q=[[Fraction(x) for x in row] for row in z['Z']]
    den=math.lcm(*(x.denominator for row in q for x in row))
    Z=np.array([[x.numerator*(den//x.denominator) for x in row] for row in q],dtype=object)
    R=np.array(proof['R'],dtype=object);n=len(Z)
    require(R.shape==Z.shape==(n,n) and np.array_equal(Z,Z.T),'Square symmetric Gram')
    require(all(R[i,i]!=0 for i in range(n)) and all(R[i,j]==0 for i in range(n) for j in range(i+1,n)),'Invertible triangular congruence')
    H=R@Z@R.T
    margin=min(int(H[i,i])-sum(abs(int(H[i,j])) for j in range(n) if i!=j) for i in range(n))
    require(margin>0 and str(margin)==proof['minimum_integer_dominance_margin'],'Strict exact dominance')
print(json.dumps({'status':'SUPPLIED_REVIEW_CONGRUENCES_RECHECK_PASS','positive_definite_grams':716,'arithmetic':'Exact integers and fractions only','optimized_python':not __debug__,'reviewer_identity_authenticated':False}))
