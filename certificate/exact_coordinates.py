#!/usr/bin/env python3
"""Exact redundant diagram coordinates for Werner31/22/13 stationarity.

No floating QR is used to define the source space or these equality maps.
A generator is the physical-site average of Herm(Q D Q), with Q the global
bosonic projector. Canonical diagram dictionaries store TOTAL orbit weight,
so a coefficient of1 means an orbit AVERAGE, not an orbit sum.
Full builds must use research-run. This script constructs no SDP or proof.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import itertools
import hashlib
import json
import math
from pathlib import Path
import resource
import signal
import time

import numpy as np
import sympy as sp
from scipy import sparse

PERMS=tuple(itertools.permutations(range(4)))
P3=tuple(itertools.permutations(range(3)))
INDEX={p:i for i,p in enumerate(PERMS)}
INV=tuple(INDEX[tuple(p.index(i) for i in range(4))] for p in PERMS)
PHYSICAL_KEEP=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23)
AUXILIARY_KEEP=(0,3,4,7,8,9,11,12,13,15,16,18,20,23)
TRIPLES=tuple(itertools.combinations_with_replacement(PHYSICAL_KEEP,3))
CATALOG=tuple(t+(a,) for t in TRIPLES for a in AUXILIARY_KEEP)
OUTPUT_TRIPLES=tuple(itertools.combinations_with_replacement(range(6),3))
OUTPUT_INDEX={t:i for i,t in enumerate(OUTPUT_TRIPLES)}
SCALE=54**4

def compose(p,q):return tuple(p[q[i]] for i in range(len(p)))
MUL=np.asarray([[INDEX[compose(p,q)] for q in PERMS] for p in PERMS],dtype=np.int64)
def cycles(p):
    seen=set();n=0
    for i in range(len(p)):
        if i not in seen:
            n+=1;j=i
            while j not in seen:seen.add(j);j=p[j]
    return n
def orbit_size(t):return math.factorial(3)//math.prod(math.factorial(t.count(v)) for v in set(t))
def canonical(word):return tuple(sorted(word[:3]))+(int(word[3]),)
def frac(x):return Fraction(int(x.p),int(x.q)) if hasattr(x,'p') else Fraction(x)
def fjson(x):return str(x.numerator) if x.denominator==1 else str(x)

def group(family):
    if family in ('31','13'):
        return tuple(INDEX[p+(3,)] for p in P3)
    if family=='22':
        return tuple(INDEX[p+q] for p in ((0,1),(1,0)) for q in ((2,3),(3,2)))
    raise ValueError(family)

def projected_generator(word,family):
    """Integer canonical-diagram numerator and common denominator2|R|².

    The physical orbit average is implicit in canonical sorting. A missing
    factor6 here would incorrectly replace the orbit average by a sum.
    """
    r=group(family);out=Counter()
    for left in r:
        for right in r:
            a=tuple(int(MUL[MUL[left,p],right]) for p in word)
            out[canonical(a)]+=1
            out[canonical(tuple(INV[p] for p in a))]+=1
    return dict(out),2*len(r)**2

def exact_local_data():
    data={}
    for d,keep in ((3,PHYSICAL_KEEP),(2,AUXILIARY_KEEP)):
        gram=sp.Matrix([[d**cycles(compose(PERMS[INV[i]],PERMS[j])) for j in range(24)] for i in range(24)])
        g=gram.extract(keep,keep);gi=g.inv();q=gi*gram.extract(keep,range(24))
        if not (gram == gram[:, list(keep)] * q): raise AssertionError('Exact consistency check failed at exact_coordinates.py:79')
        if not (q[:, list(keep)] == sp.eye(len(keep))): raise AssertionError('Exact consistency check failed at exact_coordinates.py:80')
        sparse_q=[[(i,frac(q[i,j])) for i in range(len(keep)) if q[i,j]] for j in range(24)]
        data[d]=dict(keep=keep,gram=g,inverse=gi,quotient=q,sparse_quotient=sparse_q,
                     determinant=g.det(),rank=len(keep))
    return data

def delete(p,slot,order=None):
    if order is None:order=tuple(i for i in range(4) if i!=slot)
    if sorted(order)!=[i for i in range(4) if i!=slot]:raise ValueError(order)
    out=tuple(order.index(p[slot] if p[i]==slot else p[i]) for i in order)
    embedded=list(range(4))
    for i,j in enumerate(out):embedded[order[i]]=order[j]
    return P3.index(out),INDEX[tuple(embedded)],p[slot]==slot

def auxiliary3_quotient():
    def sign(p):return (-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    q=sp.eye(5).row_join(sp.Matrix([-sign(p)//sign(P3[5]) for p in P3[:5]]))
    return [[(i,Fraction(int(q[i,j]))) for i in range(5) if q[i,j]] for j in range(6)]

class ExactMaps:
    def __init__(self):
        self.local=exact_local_data();self.aux3=auxiliary3_quotient();self.swap=INDEX[(3,1,2,0)]
        self.trace_word={d:[d**cycles(p) for p in PERMS] for d in (2,3)}
        self.objective_word={d:[d**cycles(PERMS[int(MUL[self.swap,p])]) for p in range(24)] for d in (2,3)}
        self.stiefel_local={}
        for slot in (0,3):
            dp=[];aq=[]
            for p in range(24):
                reduced,embedded,fixed=delete(PERMS[p],slot)
                dp.append((reduced,3 if fixed else 1))
                terms=Counter(dict(self.local[2]['sparse_quotient'][p]))
                for i,c in self.local[2]['sparse_quotient'][embedded]:terms[i]-=Fraction(2 if fixed else 1,2)*c
                aq.append([(i,c) for i,c in terms.items() if c])
            self.stiefel_local[slot]=(dp,aq)
        self.marginal_local={}
        for name,slot,order in [('aab',2,(0,1,3)),('abb22',1,(0,2,3)),('abb13',1,(3,0,2))]:
            dp=[];aq=[]
            for p in range(24):
                reduced,_,fixed=delete(PERMS[p],slot,order)
                dp.append((reduced,3 if fixed else 1))
                aq.append([(i,(2 if fixed else 1)*c) for i,c in self.aux3[reduced]])
            self.marginal_local[name]=(dp,aq)

    def scalar(self,word,objective=False):
        if objective:
            return math.prod(2*self.trace_word[3][p]-self.objective_word[3][p] for p in word[:3])*self.objective_word[2][word[3]]
        return math.prod(self.trace_word[3][p] for p in word[:3])*self.trace_word[2][word[3]]

    def feature(self,terms,denominator,local):
        """Output independent orbit-AVERAGE coefficients as a sparse Fraction map."""
        dp,aq=local;out=Counter()
        for word,n in terms.items():
            triple=tuple(sorted(dp[p][0] for p in word[:3]));factor=n*math.prod(dp[p][1] for p in word[:3])
            for ai,c in aq[word[3]]:out[OUTPUT_INDEX[triple]*len(self._aux_keep(local))+ai]+=factor*c
        return {i:c/denominator for i,c in out.items() if c}

    def _aux_keep(self,local):return range(14) if any(local is x for x in self.stiefel_local.values()) else range(5)

    def evaluate(self,word,family):
        terms,den=projected_generator(word,family)
        h=self.feature(terms,den,self.stiefel_local[3 if family=='13' else 0])
        result=dict(trace=Fraction(sum(n*self.scalar(w) for w,n in terms.items()),den*SCALE),
                    objective=Fraction(sum(n*self.scalar(w,True) for w,n in terms.items()),8*den*SCALE),stiefel=h)
        if family in ('31','22'):result['aab']=self.feature(terms,den,self.marginal_local['aab'])
        if family in ('22','13'):result['abb']=self.feature(terms,den,self.marginal_local['abb'+family])
        return result

    def selected_source(self,terms,den):
        """Sparse source diagram coefficients in the independent orbit-average basis."""
        out=Counter();qs=[self.local[3]['sparse_quotient']]*3+[self.local[2]['sparse_quotient']]
        for word,n in terms.items():
            for choice in itertools.product(*(q[p] for q,p in zip(qs,word))):
                key=tuple(sorted(x[0] for x in choice[:3]))+(choice[3][0],)
                out[key]+=n*math.prod(x[1] for x in choice)
        return {k:v/den for k,v in out.items() if v}

def exact_dimensions(module_path):
    archive=json.loads(module_path.read_text());characters={}
    for d in (2,3):
        case=next(c for c in archive['representations'] if c['d']==d and not c['cut'])
        characters[d]=[(m['label'],[sum(Fraction(m['basis'][p][i][i]) for i in range(m['dim'])) for p in range(24)]) for m in case['modules']]
    records={}
    for family in ('31','22','13'):
        total=0;blocks=[];r=group(family)
        for indices in itertools.combinations_with_replacement(range(len(characters[3])),3):
            stabilizer=[g for g in P3 if all(indices[g[i]]==indices[i] for i in range(3))]
            def cycle_lists(p):
                seen=set();out=[]
                for i in range(3):
                    if i not in seen:
                        cycle=[];j=i
                        while j not in seen:seen.add(j);cycle.append(j);j=p[j]
                        out.append(cycle)
                return out
            for ai,(alabel,ach) in enumerate(characters[2]):
                def character(g):
                    ans=Fraction(0)
                    for ri in r:
                        v=ach[ri]
                        for cyc in cycle_lists(g):
                            power=0
                            for _ in cyc:power=int(MUL[power,ri])
                            v*=characters[3][indices[cyc[0]]][1][power]
                        ans+=v
                    return ans/len(r)
                dimension=sum(character(g)**2+character(compose(g,g)) for g in stabilizer)/ (2*len(stabilizer))
                if not (dimension.denominator == 1): raise AssertionError('Exact consistency check failed at exact_coordinates.py:186')
                total+=int(dimension)
                if dimension:blocks.append(dict(physical_indices=indices,auxiliary_index=ai,bosonic_rank=int(character((0,1,2))),real_symmetric_invariant_dimension=int(dimension)))
        if not (total == {'31': 577, '22': 1220, '13': 577}[family]): raise AssertionError((family, total))
        records[family]=dict(exact_dimension=total,blocks=blocks)
    return records

def exact_projector_trace_dimensions(maps):
    """Independent rational rank of the idempotent on Sym^3(A3) tensor A2.

    Uses only exact local quotient coefficients, not an irrep classification.
    For K(X)=Vg X Vh (optionally after transpose), its powers still send each
    permutation diagram to one permutation diagram. Their traces are exact.
    """
    answer={}
    for family in ('31','22','13'):
        r=group(family);terms=[]
        for transpose in (False,True):
            total=Fraction(0)
            for g in r:
                for h in r:
                    traces={}
                    for d in (3,2):
                        keep=maps.local[d]['keep'];q=maps.local[d]['quotient'];power_traces=[]
                        labels=list(keep)
                        for power in range(1,4):
                            labels=[int(MUL[MUL[g,INV[p] if transpose else p],h]) for p in labels]
                            power_traces.append(sum(frac(q[i,p]) for i,p in enumerate(labels)))
                        traces[d]=power_traces
                    t1,t2,t3=traces[3]
                    total+=(t1**3+3*t1*t2+2*t3)*traces[2][0]/6
            terms.append(total)
        dimension=sum(terms)/(2*len(r)**2)
        if not (dimension == {'31': 577, '22': 1220, '13': 577}[family]): raise AssertionError('Exact consistency check failed at exact_coordinates.py:219')
        answer[family]=dict(exact_dimension=int(dimension),sandwich_trace_sum=fjson(terms[0]),
            transpose_sandwich_trace_sum=fjson(terms[1]),projector_denominator=2*len(r)**2)
    return answer

def pack_fraction_matrix(columns,rows):
    den=math.lcm(*(x.denominator for col in columns for x in col.values())) if any(columns) else 1
    rr=[];cc=[];vv=[]
    for j,col in enumerate(columns):
        for i,x in col.items():rr.append(i);cc.append(j);vv.append(x.numerator*(den//x.denominator))
    if not (max(map(abs, vv), default=0) <= np.iinfo(np.int64).max): raise AssertionError('Exact consistency check failed at exact_coordinates.py:229')
    m=sparse.csc_matrix((np.asarray(vv,dtype=np.int64),(rr,cc)),shape=(rows,len(columns)))
    return m,den

def modular_rank(matrix,prime=1000003):
    """Exact finite-field rank, only a LOWER BOUND on rational rank."""
    basis={}
    for j in range(matrix.shape[1]):
        v=np.asarray(matrix[:,j].toarray(),dtype=np.int64).ravel()%prime
        for pivot,b in basis.items():
            if v[pivot]:v=(v-v[pivot]*b)%prime
        nz=np.flatnonzero(v)
        if len(nz):
            p=int(nz[0]);v=(v*pow(int(v[p]),prime-2,prime))%prime;basis[p]=v
    return len(basis)
