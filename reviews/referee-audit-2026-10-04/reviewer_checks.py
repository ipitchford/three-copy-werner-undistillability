#!/usr/bin/env python3
"""Alternative exact arithmetic replay. Author coordinate/foundation definitions retained.
The PSD proof uses integer congruence and strict diagonal dominance, not floating
Cholesky as an acceptance test. Float arithmetic ONLY proposes an integer matrix.
Run from a clean output directory; no author replay receipts are trusted.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
os.environ['SYMPY_GROUND_TYPES']='python'
import sys,math,json,hashlib,time,types,traceback
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
import numpy as np
import sympy as sp
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
ROOT=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/werner_review/werner-problem6-proof').resolve()
OUT=Path(sys.argv[2] if len(sys.argv)>2 else '/mnt/data/werner_referee_audit_2026-10-04').resolve()
OUT.mkdir(exist_ok=True,parents=True)
C=ROOT/'certificate';sys.path.insert(0,str(C))
START=time.monotonic()
def log(*x):print(round(time.monotonic()-START,3),*x,flush=True)
def write(name,x):(OUT/name).write_text(json.dumps(x,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(c,m):
 if not c:raise AssertionError(m)
class NMod:
 """Exact prime-field determinant; multiplication <=(p-1)^2 fits int64 here."""
 def __init__(self,rows,p):
  self.p=int(p);require(self.p<=1000003,'Modulus integer safety');self.a=np.asarray(rows,dtype=np.int64)%self.p
 def det(self):
  a=self.a.copy();p=self.p;n=len(a);require(a.shape==(n,n),'Square minor');d=1
  for k in range(n):
   piv=np.flatnonzero(a[k:,k])
   if not len(piv):return 0
   j=k+int(piv[0])
   if j!=k:a[[k,j]]=a[[j,k]];d=-d
   v=int(a[k,k]);d=d*v%p
   if k+1<n:
    f=a[k+1:,k]*pow(v,-1,p)%p
    a[k+1:,k+1:]=(a[k+1:,k+1:]-f[:,None]*a[k,k+1:][None,:])%p
  return d%p
shim=types.ModuleType('flint');shim.nmod_mat=NMod
sys.modules['flint']=shim
from foundation import verify_foundation
from basis_values import basis_values
from exact_coordinates import PERMS,compose

def clear(rows):
 a=np.asarray(rows,dtype=object);vals=[Fraction(x) for x in a.ravel()];den=math.lcm(*(v.denominator for v in vals))
 return np.asarray([v.numerator*(den//v.denominator) for v in vals],dtype=object).reshape(a.shape),den

def norm(a,d):
 g=math.gcd(d,*(int(v) for v in np.asarray(a).ravel()))
 if g>1:return a//g,d//g
 return a,d

def psd_exact(rows):
 z,zd=clear(rows);n=len(z);require(z.shape==(n,n) and np.array_equal(z,z.T),'Symmetric Z')
 if n==0:return {'order':0,'PD':True,'R':[]},z,zd
 zz=np.asarray([[float(Fraction(v)) for v in row] for row in rows])
 ll=np.linalg.cholesky(zz);inv=np.linalg.inv(ll)
 for scale in (10**5,10**8,10**11,10**14):
  r=np.asarray([[int(round(float(v)*scale)) for v in row] for row in inv],dtype=object)
  require(all(r[i,j]==0 for i in range(n) for j in range(i+1,n)),'Lower triangular congruence')
  if any(r[i,i]==0 for i in range(n)):continue
  h=(r@z)@r.T
  margins=[int(h[i,i])-sum(abs(int(h[i,j])) for j in range(n) if j!=i) for i in range(n)]
  if min(margins)>0:
   return {'order':n,'PD':True,'method':'exact integer invertible congruence to positive strictly diagonally dominant matrix',
    'R':r.tolist(),'minimum_integer_dominance_margin':str(min(margins)),'Z_denominator_bits':zd.bit_length()},z,zd
 raise AssertionError('No diagonal-dominance certificate found; no assertion of indefiniteness')

ARCH=None
@lru_cache(maxsize=None)
def local_factor(ci,mi,label):
 module=ARCH['representations'][ci]['modules'][mi]
 return clear([[x for row in module['basis'][p] for x in row] for p in label])
@lru_cache(maxsize=None)
def gram(refs):
 g=np.ones((1,1),dtype=object);den=1
 for ci,mi in refs:
  a,d=clear(ARCH['representations'][ci]['modules'][mi]['gram']);g=np.kron(g,a);den*=d
 return norm(g,den)
def table(weight,den,refs,labels):
 dims=[ARCH['representations'][ci]['modules'][mi]['dim'] for ci,mi in refs]
 value=weight.T.reshape(tuple(dims+dims)).transpose(0,4,1,5,2,6,3,7).reshape(tuple(d*d for d in dims))
 factors=[]
 for ci,mi in refs:
  f,d=local_factor(ci,mi,tuple(labels));factors.append(f);den*=d
 for axis in sorted(range(4),key=lambda i:dims[i],reverse=True):
  value=np.moveaxis(np.tensordot(factors[axis],value,axes=(1,axis)),0,axis)
 return value,den

def run():
 global ARCH
 # Content hashes protect the exact submission being reviewed, not correctness.
 m=json.loads((ROOT/'MANIFEST.json').read_text())
 for name,entry in m['files'].items():
  p=ROOT/name;require(p.resolve().is_relative_to(ROOT),'Safe path')
  require(p.stat().st_size==entry['bytes'] and sha(p)==entry['sha256'],'Top-level hash '+name)
 cm=json.loads((C/'manifest.json').read_text())
 for name,digest in cm['files'].items():require(sha(C/name)==digest,'Certificate hash '+name)
 write('hash_check.json',{'top_level_files':len(m['files']),'certificate_files':len(cm['files']),'all_hashes_match':True,
 'certificate_manifest_sha256':sha(C/'manifest.json'),'input_zip_sha256':sha(ROOT.parent.parent/'werner-problem6-proof.zip')})
 log('All hashes match',len(m['files']),len(cm['files']))
 foundation=verify_foundation(C)
 foundation['reviewer_backend']='Author foundation/coordinate logic with independent numpy integer prime-field determinant; SymPy exact rationals; no FLINT.'
 write('foundation.json',foundation);log('Foundation pass')
 ARCH=json.loads((C/'data/local_modules.json').read_text())
 cert=json.loads((C/'certificate.json').read_text());recs=cert['retained_multipliers']
 with np.load(C/'data/basis_data.npz',allow_pickle=False) as data:
  indices={f:data['basis'+f+'_indices'] for f in ['31','22','13']}
  B=data['row_basis_integer'].astype(object);obj=data['global_objective_numerator'].astype(object);od=int(data['global_objective_denominator']);eqrows=data['equality_basis_rows'].tolist()
 widths={'31':577,'22':1220,'13':577};offsets={'31':0,'22':577,'13':1797}
 tau=(0,1,3,2);cross=(2,3,1,0)
 diagonal=[PERMS.index(compose(tau,compose(p,tau))) for p in PERMS]
 crossmap=[PERMS.index(compose(p,cross)) for p in PERMS]
 total=np.zeros(2374,dtype=object);td=1;proofs=[];compared=0
 for i,rec in enumerate(recs):
  block=json.loads((C/rec['path']).read_text());require(block['id']==rec['id'],'Block identity')
  refs=tuple((int(r['case']),int(r['module'])) for r in block['modules']);dims=[]
  for ref,(ci,mi),d in zip(block['modules'],refs,[3,3,3,2]):
   case=ARCH['representations'][ci];module=case['modules'][mi]
   require(case['d']==d==ref['d'] and case['cut']==block['cut'],'Cut/dimension binding')
   require(module['dim']==ref['dim'] and module['label']==ref['label'],'Module binding');dims.append(module['dim'])
  proof,z,zd=psd_exact(block['Z']);proof['id']=rec['id'];proofs.append(proof)
  N,nd=clear(block['N']);require(list(N.shape)==block['N_shape'] and list(z.shape)==block['Z_shape'] and N.shape[1]==len(z),'Gram shapes')
  n=math.prod(dims);g,gd=gram(refs)
  if block['kind']=='source':
   require(N.shape[0]==n,'Source shape');ops=[(block['family'],N,N,list(range(24)),1)]
  else:
   require(block['kind']=='joint11' and block['cut']==[3] and N.shape[0]==2*n,'Coupled shape')
   top=N[:n];bottom=N[n:];ops=[('31',top,top,diagonal,1),('22',bottom,top,crossmap,2),('13',bottom,bottom,diagonal,1)]
  # Compute from N,Z, modules and basis first; saved author functionals only compare after.
  with np.load(ROOT/'replay'/('block_'+str(i).zfill(4)+'.npz'),allow_pickle=False) as old:
   for f,left,right,mapping,factor in ops:
    w,wd=norm((factor*(left@z@right.T))@g,nd*nd*zd*gd)
    tab,d=table(w,wd,refs,mapping)
    values=basis_values(tab,d,f,indices[f]);v,vd=clear(values);v,vd=norm(v,vd)
    del tab,w
    require(len(v)==widths[f],'All basis entries evaluated')
    ov=np.asarray([int(x) for x in old['functional'+f+'_numerator']],dtype=object);ovd=int(old['functional'+f+'_denominator'])
    require(np.array_equal(v*ovd,ov*vd),'Author functional differs '+str(i)+' '+f);compared+=1
    common=math.lcm(td,vd);total*=common//td;start=offsets[f];total[start:start+len(v)]+=v*(common//vd);td=common
  if (i+1)%20==0 or i==0:log('Blocks',i+1,'/',len(recs),'exact PSD and full functionals pass')
 equality=json.loads((C/'data/equality_multiplier.json').read_text());require(equality['equality_basis_rows']==eqrows,'Selected equality rows')
 eta,ed=clear(equality['eta']);eq=B.T@eta
 common=math.lcm(td,ed,od);residual=obj*(common//od)-total*(common//td)-eq*(common//ed)
 require(not any(residual),'Global exact residual')
 write('psd_congruence_certificates.json',{'proofs':proofs})
 summary={'status':'REVIEWER_ALTERNATIVE_EXACT_ARITHMETIC_PASS','all_hashes_match':True,'foundation_recomputed':True,
 'multipliers_recomputed':len(recs),'strictly_positive_definite_small_Z':len(proofs),'source_multiplier_count':sum(x['kind']=='source' for x in recs),
 'coupled_multiplier_count':sum(x['kind']=='joint11' for x in recs),'functional_vectors_compared':compared,
 'exact_global_equations':2374,'all_integer_residuals_zero':True,'common_denominator_bits':common.bit_length(),
 'elapsed_seconds':time.monotonic()-START,'basis_and_foundation_logic':'Author definitions retained; not a fully independent second mathematical implementation',
 'independent_backends':['Prime-field determinant','Rational integer tensor contraction','PSD congruence and diagonal-dominance proof','Global rational summation'],
 'native_author_FLINT_verifier_run':False,'reason_native_not_run':'python-flint unavailable; installation/network download unsuccessful',
 'saved_author_receipts_trusted':False}
 write('summary.json',summary);log(summary)
try:run()
except Exception as e:
 write('failure.json',{'exception':repr(e),'traceback':traceback.format_exc(),'elapsed_seconds':time.monotonic()-START});traceback.print_exc();sys.exit(1)
