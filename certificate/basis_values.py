"""Cached exact projected-generator evaluation for integer24^4 trace tables.

Drop-in replacement for werner_exact_joint_functional_table.basis_values.
Only integer contraction is vectorized. All sums use Python object integers;
there is no fixed-width multiplication or floating conversion in the value.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import itertools
from numbers import Integral

import numpy as np
from exact_coordinates import CATALOG,projected_generator


@dataclass(frozen=True)
class CompiledBasis:
    family: str
    generator_indices: tuple
    table_indices: np.ndarray
    coefficients: np.ndarray
    row_starts: np.ndarray
    row_denominators: tuple

    def values(self,table,denominator):
        if isinstance(denominator,bool) or not isinstance(denominator,Integral) or denominator<=0:
            raise TypeError('Trace-table denominator must be a positive exact integer')
        array=np.asarray(table)
        if array.shape!=(24,24,24,24):raise ValueError('Expected24^4 integer trace table')
        if array.dtype.kind not in 'iuO':raise TypeError('Trace table must contain integers')
        if not self.generator_indices:return []
        # Force object BEFORE multiplication, even when the incoming table is
        # int64. Its exact entries and coefficients can have arbitrarily large
        # products or sums. Integral result checks reject object-float tables.
        selected=array.ravel()[self.table_indices].astype(object,copy=False)
        products=selected*self.coefficients
        sums=np.add.reduceat(products,self.row_starts)
        if any(isinstance(x,bool) or not isinstance(x,Integral) for x in sums):
            raise TypeError('Noninteger values in trace-table contraction')
        return [Fraction(int(n),int(denominator)*d) for n,d in zip(sums,self.row_denominators)]


@lru_cache(maxsize=32)
def _compile(family,indices):
    if family not in ('31','22','13'):raise ValueError('Unknown source family')
    if any(i<0 or i>=len(CATALOG) for i in indices):raise IndexError('Generator index outside catalog')
    flat=[];coefficients=[];starts=[];denominators=[]
    for index in indices:
        starts.append(len(flat));terms,den=projected_generator(CATALOG[index],family)
        for word,n in terms.items():
            permutations=set(itertools.permutations(word[:3]));factor=int(n)*(6//len(permutations))
            for p in permutations:
                flat.append(((p[0]*24+p[1])*24+p[2])*24+word[3])
                coefficients.append(factor)
        if len(flat)==starts[-1]:raise ValueError('Unexpected empty projected-generator row')
        denominators.append(6*int(den))
    if not (all((d == (432 if family in ('31', '13') else 192) for d in denominators))): raise AssertionError('Exact consistency check failed at basis_values.py:58')
    ix=np.asarray(flat,dtype=np.int64);cs=np.asarray(coefficients,dtype=object)
    rs=np.asarray(starts,dtype=np.int64)
    for a in (ix,cs,rs):a.flags.writeable=False
    return CompiledBasis(family,indices,ix,cs,rs,tuple(denominators))


def compile_basis(family,indices):
    converted=[]
    for i in indices:
        if isinstance(i,bool) or not isinstance(i,Integral):raise TypeError('Generator indices must be exact integers')
        converted.append(int(i))
    return _compile(str(family),tuple(converted))


def basis_values(table,denominator,family,indices):
    return compile_basis(family,indices).values(table,denominator)
