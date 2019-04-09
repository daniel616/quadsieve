
import math
import numpy as np

import quadsieve as qs

#This code factors n from start to finish.
def solve(n):
    #breakpoint()
    fact_mat,factors,row_labels,row_labels_unsquared=qs.quad_sieve(n)

    row_labels=np.array(row_labels,dtype="int64")
    row_labels_unsquared=np.array(row_labels_unsquared,dtype="int64")
    factors=np.array(factors,dtype="int64")
    to_square=np.array(row_labels_unsquared,dtype="int64")
    print("fact_mat.shape",fact_mat.shape)
    print("row_labels shape", row_labels.shape)
    print("factors shape", factors.shape)
    print("row_labels_unsquared shape", to_square.shape)

    parsed = np.mod(fact_mat, 2 * np.ones(fact_mat.shape)).astype('bool')
    transform, reduced = row_echelon(parsed)

    null_row_inds=[]
    for idx, row in enumerate(reduced):
        if check_null(row):
            null_row_inds.append(idx)

    for x in null_row_inds:
        null_row = transform[x]
        val=get_factor(null_row,to_square,fact_mat,factors,n)
        if val>1 and val<n:
            return val,int(n/val)


    return None

#Takes an array of factors and exponents and as input.
#Returns product of factors raised to exponents mod n.
def factors_to_num(factors,exps, n):
    ret = 1
    for idx,_ in enumerate(factors):
        power=pow(int(factors[idx]),int(exps[idx]),n)
        ret= (ret * power)%n
    return ret

#Returns true iff all elements of row are 0.
def check_null(row):
    for x in row:
        if x!=0: return False
    return True

#
def get_factor(null_row,to_square,fact_mat,factors,n):
    print(null_row.shape)
    print(to_square.shape)
    assert null_row.shape==to_square.shape

    exps=np.matmul(null_row,fact_mat)

    for idx, x in np.ndenumerate(exps):
        assert x%2==0

    exps=exps/2

    sqrt_val=factors_to_num(factors,exps,n)

    orig_val=1
    for idx, x in np.ndenumerate(null_row):
        if x==1:
            orig_val=(orig_val*to_square[idx])%n

    print("orig_val:",orig_val,"sqrt_val:",sqrt_val)
    orig_val,sqrt_val=int(orig_val),int(sqrt_val)
    assert (orig_val**2)%n==(sqrt_val**2)%n
    ret= math.gcd(orig_val-sqrt_val,n)
    print(ret)
    return ret


#Input: 2d np boolean array. Treated as Z/2Z. Returns reduced row echelon form of array, and matrix of row
#transformations
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
#debugging purposes only.
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

if __name__ == '__main__':
    #print(solve(101*61))
    #print(solve(10_172_605_169))
    print(solve(1000000007*1000000009))
    #print(solve(16921456439215439701))
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

