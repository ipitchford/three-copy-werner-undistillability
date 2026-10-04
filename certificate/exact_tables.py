"""Portable exact PSD and tensor-table arithmetic. No floating arithmetic."""
from fractions import Fraction
from functools import lru_cache
import hashlib
import math
import numpy as np
from flint import fmpq_mat
from exact_coordinates import PERMS,compose
from basis_values import basis_values

def require(value,message):
    if not value:raise AssertionError(message)

def qmatrix(rows):
    return fmpq_mat(rows)

def qblock(matrix,rows,columns):
    return fmpq_mat([[matrix[i,j] for j in columns] for i in rows])

def qkron(matrices):
    result=fmpq_mat([[1]])
    for matrix in matrices:
        a,b=result.nrows(),matrix.nrows();out=fmpq_mat(a*b,a*b)
        for i in range(a):
            for j in range(a):
                value=result[i,j]
                if not value:continue
                for k in range(b):
                    for l in range(b):
                        if matrix[k,l]:out[i*b+k,j*b+l]=value*matrix[k,l]
        result=out
    return result

def charpoly_psd(matrix):
    """For a symmetric rational matrix, alternating signs exactly test PSD.

    All eigenvalues are real. Alternating nonnegative characteristic
    coefficients imply (-1)^n p(-t)>0 for t>0 and hence no negative root.
    The converse follows from the elementary symmetric functions.
    """
    n=matrix.nrows();require(n==matrix.ncols(),'Square Gram')
    require(matrix==matrix.transpose(),'Exactly symmetric Gram')
    if not n:return dict(PSD=True,PD=True,rank=0,order=0)
    integer,den=matrix.numer_denom();require(den>0,'Positive denominator')
    polynomial=integer.charpoly();coefficients=[int(polynomial[i]) for i in range(n+1)]
    require(coefficients[-1]==1,'Monic characteristic polynomial')
    require(all(((-1)**k)*coefficients[n-k]>=0 for k in range(n+1)),'Negative Gram eigenvalue')
    nullity=next(i for i,c in enumerate(coefficients) if c)
    return dict(PSD=True,PD=nullity==0,rank=n-nullity,order=n,
        method='symmetric integer characteristic polynomial signs',
        clearing_denominator_bits=int(den).bit_length(),maximum_coefficient_bits=max(abs(c).bit_length() for c in coefficients))

def integer_array(matrix):
    integer,den=matrix.numer_denom();require(den>0,'Positive matrix clearing denominator')
    return np.asarray([[int(integer[i,j]) for j in range(integer.ncols())] for i in range(integer.nrows())],dtype=object),int(den)

def local_factor(module,label_map):
    key=tuple(map(int,label_map));cache=module.setdefault('_factors',{})
    if key not in cache:
        rows=[[x for row in module['basis'][p] for x in row] for p in key]
        cache[key]=integer_array(qmatrix(rows))
    return cache[key]

def trace_table(weight,modules,label_map):
    """T[p]/den = Tr(weight * tensor_i rho_i(label_map[p_i]))."""
    dims=[int(m['dim']) for m in modules]
    integers,den=integer_array(weight.transpose())
    require(integers.shape==(math.prod(dims),)*2,'Tensor weight shape')
    value=integers.reshape(tuple(dims+dims)).transpose((0,4,1,5,2,6,3,7)).reshape(tuple(d*d for d in dims))
    factors=[]
    for module in modules:
        factor,fd=local_factor(module,label_map);factors.append(factor);den*=fd
    for axis in sorted(range(4),key=lambda i:dims[i],reverse=True):
        value=np.moveaxis(np.tensordot(factors[axis],value,axes=(1,axis)),0,axis)
    require(value.shape==(24,)*4,'Trace table shape')
    return value,den

def permutation_maps():
    tau=(0,1,3,2);cross=(2,3,1,0)
    diagonal=[PERMS.index(compose(tau,compose(p,tau))) for p in PERMS]
    off_diagonal=[PERMS.index(compose(p,cross)) for p in PERMS]
    return diagonal,off_diagonal

def block_functionals(block,modules,basis_indices,deadline=lambda:None):
    """Replay final N,Z directly. Returns PSD proof and exact family vectors."""
    n=math.prod(int(m['dim']) for m in modules);N=qmatrix(block['N']);Z=qmatrix(block['Z'])
    require([N.nrows(),N.ncols()]==block['N_shape'],'N declared shape')
    require([Z.nrows(),Z.ncols()]==block['Z_shape'],'Z declared shape')
    require(Z.nrows()==N.ncols(),'Gram factor dimensions');proof=charpoly_psd(Z);deadline()
    G=qkron([qmatrix(m['gram']) for m in modules]);vectors={}
    if block['kind']=='source':
        require(N.nrows()==n and block['family'] in ('31','22','13'),'Source dimension/family')
        operations=[(block['family'],N,N,list(range(24)),1)]
    else:
        require(block['kind']=='joint11' and block['cut']==[3] and N.nrows()==2*n,'Joint convention/dimension')
        top=qblock(N,range(n),range(N.ncols()));bottom=qblock(N,range(n,2*n),range(N.ncols()))
        diagonal,cross=permutation_maps()
        # Do NOT symmetrize the lower-left cross multiplier.
        operations=[('31',top,top,diagonal,1),('22',bottom,top,cross,2),('13',bottom,bottom,diagonal,1)]
    for family,left,right,mapping,factor in operations:
        weight=factor*(left*Z*right.transpose())*G
        table,den=trace_table(weight,modules,mapping);deadline()
        vectors[family]=basis_values(table,den,family,basis_indices[family])
        del table,weight;deadline()
    return proof,vectors

def clear_vector(values):
    den=math.lcm(*(v.denominator for v in values))
    numbers=[v.numerator*(den//v.denominator) for v in values]
    return numbers,den
