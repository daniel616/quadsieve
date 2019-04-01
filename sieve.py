
from math import *
import numpy as np

import quadsieve as qs


def solve(n):
    fact_mat,factors,to_square=qs.quad_sieve(n)
    print(factors)
    print(to_square)
    print(fact_mat)
    parsed = np.mod(fact_mat, 2 * np.ones(fact_mat.shape)).astype('bool')

    transform, reduced = row_echelon(parsed)

    #print(transform)
    #print(reduced)

    null_row = transform[transform.shape[0] - 1]

    s_exps = np.matmul(null_row, fact_mat)


def get_factor(null_row,fact_mat,factors,to_square,n):
    pass

#Input: 2d np boolean array. Treated as Z/2Z. Returns reduced row echelon form of array, and matrix of transformations
#used to obtain this form.
def row_echelon(mat):
    orig=mat
    mat=np.array(mat,copy=True)

    mir=np.identity(mat.shape[0],dtype="bool")
    nrows,ncols=mat.shape[0],mat.shape[1]


    r_base=0
    for j in range(ncols):
        if r_base==nrows: break
        for i in range(r_base,nrows):
            if mat[i][j]==True:
                switch_row(mat,i,r_base,mirror=mir)
                clear_pivot(mat,r_base,j,mirror=mir)
                r_base+=1
                break

    #print(matmul2(mir,orig),"product\n")
    #print(mat,"final\n")
    assert np.array_equal(matmul2(mir,orig),mat)
    return mir, mat

#Multiplies matrices over Z/2Z. Regular matmul doesn't work because + is not xor.
def matmul2(a,b):
    c_rows=[]
    for i in range(a.shape[0]):
        row_i=np.zeros([b.shape[1]])
        for j in range(a.shape[1]):
            if a[i][j]==1:
                row_i=np.logical_xor(row_i,b[j])
        c_rows.append(row_i)
    return np.array(c_rows)

#Selects a pivot in a matrix over 2Z, and subtracts pivot
# row from rows below to make sure
#that there are no nonzero elements below the pivot.
def clear_pivot(mat,i,j,mirror=None):
    for r in range(i+1,mat.shape[0]):
        if mat[r][j]==True:
            mat[r]=np.logical_xor(mat[r],mat[i])
            if mirror is not None:
                assert mirror.shape[0]==mat.shape[0]
                mirror[r]=np.logical_xor(mirror[r],mirror[i])


def switch_row(mat, r1,r2,mirror=None):
    tmp=np.array(mat[r1],copy=True)
    mat[r1]=mat[r2]
    mat[r2]=tmp

    if mirror is not None:
        assert mirror.shape[0]==mat.shape[0]
        switch_row(mirror,r1,r2)

def gcd(a,b):
    if a==0: return b
    if b==0: return a

    if a>=b: return gcd(a%b,b)
    else: return gcd(a,b%a)

if __name__ == '__main__':
    solve(101*67)
    '''
    mat1=np.array([[1,0,1,1],
                   [0,1,0,1],
                   [1,0,0,0],
                   [1,1,0,0],
                   [1,0,1,0]],dtype="bool")
    print(mat1,"\n")
    mir,red=row_echelon(mat1)
    print(red,"\n")
    print(mir)'''
