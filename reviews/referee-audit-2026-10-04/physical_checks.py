#!/usr/bin/env python3
"""Independent small exact Gaussian-integer tests of physical conventions.
Finite tests validate conventions, not the universal rank-two inequality.
"""
from pathlib import Path
from fractions import Fraction
import itertools,json
import numpy as np
OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261004)
def conj(x):return (x[0],-x[1])
def mult(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def prod(xs):
 v=(1,0)
 for x in xs:v=mult(v,x)
 return v
def trace(C,S):
 t=C.reshape((3,)*6);n=3
 for i in sorted(S,reverse=True):t=np.trace(t,axis1=i,axis2=i+n);n-=1
 return t
SUBSETS=[tuple(i for i in range(3) if s>>i&1) for s in range(8)]
def norm2(re,im):return sum(int(v)**2 for v in np.asarray(re).ravel())+sum(int(v)**2 for v in np.asarray(im).ravel())
def q3(re,im):
 return sum((Fraction(-1,2)**len(s))*norm2(trace(re,s),trace(im,s)) for s in SUBSETS)
omega=np.zeros(9,dtype=object);omega[[0,4,8]]=1
w=2*np.eye(9,dtype=object)-np.outer(omega,omega);W=np.kron(np.kron(w,w),w)
zero=np.zeros((3,3),dtype=object);zero[0,0]=1
D=np.diag([1,1,0]).astype(object);E01=np.zeros((3,3),dtype=object);E01[0,1]=1
Z1=np.kron(np.kron(D,zero),zero);Z2=np.kron(np.kron(D,E01),zero);R3=np.kron(np.kron(np.eye(3,dtype=object),zero),zero)
checks={'normal_rank_two_zero':str(q3(Z1,np.zeros_like(Z1))),
'nonnormal_rank_two_zero':str(q3(Z2,np.zeros_like(Z2))),
'rank_three_negative':str(q3(R3,np.zeros_like(R3)))}
assert checks=={'normal_rank_two_zero':'0','nonnormal_rank_two_zero':'0','rank_three_negative':'-3/8'}
swap_tests=0;direct_tests=0;non_normal=0;coupled_tests=0
for case in range(24):
 r=1 if case<8 else 2
 ur=rng.integers(-2,3,(27,r)).astype(object);ui=rng.integers(-2,3,(27,r)).astype(object)
 vr=rng.integers(-2,3,(27,r)).astype(object);vi=rng.integers(-2,3,(27,r)).astype(object)
 cr=ur@vr.T+ui@vi.T;ci=ui@vr.T-ur@vi.T
 ar=np.zeros((27,2),dtype=object);ai=ar.copy();br=ar.copy();bi=ar.copy()
 ar[:,:r]=ur;ai[:,:r]=ui;br[:,:r]=vr;bi[:,:r]=vi
 # Check q3 against direct W^{tensor3}; tensor slots interleave A_i,B_i.
 pr=cr.reshape((3,)*6).transpose(0,3,1,4,2,5).ravel();pi=ci.reshape((3,)*6).transpose(0,3,1,4,2,5).ravel()
 direct=Fraction(int(pr@(W@pr)+pi@(W@pi)),8)
 assert direct==q3(cr,ci);assert direct>=0;direct_tests+=1
 # Real part of CC* minus C*C, sufficient here to detect nonnormal examples.
 non_normal+=not np.array_equal(cr@cr.T+ci@ci.T,cr.T@cr+ci.T@ci)
 tr=np.outer(ar.ravel(),br.ravel())-np.outer(ai.ravel(),bi.ravel())
 ti=np.outer(ar.ravel(),bi.ravel())+np.outer(ai.ravel(),br.ravel())
 tr=tr.reshape(3,3,3,2,3,3,3,2);ti=ti.reshape(tr.shape)
 for s in SUBSETS:
  axes=list(range(8))
  for j in (*s,3):axes[j],axes[j+4]=axes[j+4],axes[j]
  value=int(np.sum(tr*tr.transpose(axes)+ti*ti.transpose(axes)))
  assert value==norm2(trace(cr,s),trace(ci,s));swap_tests+=1
 # Literal nonsymmetrized coupled realignment: 10 random entries per case.
 a=[(int(x),int(y)) for x,y in zip(ar.ravel(),ai.ravel())];b=[(int(x),int(y)) for x,y in zip(br.ravel(),bi.ravel())]
 for _ in range(10):
  i=rng.integers(0,54,4).tolist();j=rng.integers(0,54,4).tolist()
  row=(i[0],i[1],i[2],j[3]);col=(i[3],j[2],j[0],j[1]);fs=(a,a,b,b)
  lhs=prod([term for f,x,y in zip(fs,row,col) for term in (f[x],conj(f[y]))])
  u=prod([a[i[0]],a[i[1]],b[i[2]],conj(a[i[3]])]);v=prod([b[j[0]],b[j[1]],a[j[2]],conj(b[j[3]])])
  assert lhs==mult(u,conj(v));coupled_tests+=1
checks.update(status='EXACT_PHYSICAL_CONVENTION_TESTS_PASS',arithmetic='Python arbitrary-precision Gaussian integers and fractions',seed=20261004,
 arbitrary_complex_rank_at_most_two_cases=24,direct_W_tensor_checks=direct_tests,auxiliary_swap_partial_trace_checks=swap_tests,
 literal_nonsymmetric_cross_map_checks=coupled_tests,detected_nonnormal_cases=int(non_normal),
 limitation='Finite calibration tests; not a universal proof and not an independent coordinate-basis implementation')
(OUT/'physical_checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
