
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import math as m
import numpy as num

def gen_primes(n):
    #generates primes up to n using sieve of eratosthenes
    
    rit = list(range(2, n + 1))
    primes = []
    
    for i in range(len(rit)):
        #any primes won't be factored so they'll = 1 in the list
        if (rit[i] != 1):
            #each composite number gets divided by its prime factors
            for j in range(i + rit[i], len(rit), rit[i]):
                while(rit[j] % rit[i] == 0):
                    rit[j] //= rit[i]
    
    #any number != 1 is prime
    for val in rit:
        if (val != 1):
            primes.append(val)
    
    return primes

def quad_sieve(n):
    s = isqrt(n) + 1 #low x bound
    B = gen_primes(m.ceil((1/2.2)*m.exp(m.sqrt(m.log(n)*m.log(m.log(n)))))) #formula for most efficient B
    upx = m.ceil(m.exp((2/2.2)*m.sqrt(m.log(n)*m.log(m.log(n))))) #formula for upper x bound
    xrng = list(range(s, s + upx))
    vals = [x**2 - n for x in xrng] #generates list of x^2 - n values
    valsc = vals.copy() #copy of x^2 - n values, once you have indices of B-smooth numbers you can recover the number
    indices = [] #indices of B-smooth numbers in vals, valsc
    new_B = [] #discard primes in B that aren't squares mod n
    factor_found = 0 #b-smooth numbers found
    print('hi')
    #discard non squares p mod n
    for p in B:
        if p == 2:
            new_B.append(p)
        else:
            if (legendre_symbol(n, p) == 1):
                new_B.append(p)
    
    #returned matrix will have 1 more factor than there are primes in new_B to gurantee a linear dependency
    factor_matrix = num.zeros((len(new_B) + 1, len(new_B)), dtype=int)
    
    #performs sieve of eratosthenes on x^2 - n sequence for primes in new_B
    for p in new_B:
        #special case when p = 2, find first x such that x^2 - n divisible by 2 then every other x after will be divisible by 2
        if p == 2:
            i = 0;
            
            while vals[i] % 2 != 0:
                i += 1
            
            for j in range(i, len(vals), 2):
                #account for prime powers
                while(vals[j] % 2 == 0):
                    vals[j] //= 2
                
                if(vals[j] == 1):
                    factor_found += 1
                    #break if enough b-smooth numbers are found
                    if (factor_found == len(new_B) + 1):
                        break
        else:
            #solve x^2 = n mod p, two values a1 and a2
            avals = modular_sqrt(n, p)
            
            for a in avals:
                #find first x = a1, a2 mod p
                for x in range(s, s + upx):
                    if (x % p == a):
                        i = x - s
                        break
                #x+p, x+2p, x+3p, ... = a1, a2 mod p also so x^2 - n is divisible by p
                for j in range(i, len(vals), p):
                    #divide by p at these x values, account for prime powers
                    while(vals[j] % p == 0):
                        vals[j] //= p
                    
                    #every b-smooth number will equal 1
                    if(vals[j] == 1):
                        factor_found += 1
                        if (factor_found == len(new_B) + 1):
                            break
                
                if (factor_found == len(new_B) + 1):
                    break
        
        if (factor_found == len(new_B) + 1):
            break;
    
    #find indices in vals where vals[i] == 1, so vals[i] is b-smooth
    for i in range(len(vals)):
        if (vals[i] == 1):
            indices.append(i)


    # print(new_B)
    countr = 0
    row_titles = []
    for index in indices:
        factors = {key : 0 for key in new_B}
        row_title = valsc[index]
        if row_title != 0:
            print("Creating a new row: %d" % row_title)
            row_titles.append(row_title)
            new_row = [int(factors[x]) for x in sorted(factor(row_title,new_B,factors))]
            factor_matrix[countr] = new_row
            countr+=1
    # print(factor_matrix)
    # print(len(new_B))
    # print(factor_found)

    return (factor_matrix,new_B,row_titles,[xrng[indexex] for indexex in indices])

def factor(n,new_B,factors):
    if factors==-1:
        return -1
    #print("going to factor %d" % n)
    B = new_B[-1]
    toFactor = n
    if toFactor in new_B:
        factors[toFactor] += 1
        #print("Found factor! Switching %d" % toFactor)
        return factors

    nontrivial_factor = pollard_rho(toFactor)
    if nontrivial_factor > 0:
        #print(nontrivial_factor)
        if nontrivial_factor in new_B:
            factors[nontrivial_factor] += 1
            #print("Found factor! Switching %d" % nontrivial_factor)
            factors = factor(n//nontrivial_factor,new_B,factors)
            return factors
        else:
            factors = factor(nontrivial_factor,new_B,factors)
            factors = factor(n//nontrivial_factor,new_B,factors)

    else:  
        #print("pollard_rho couldn't factor %d" % n)
        return -1
    return factors

def pollard_rho(n):
    for c in range(1,6):
        if n % 2 == 0:
            return n//2        
        print("finding nontrival factor of %d" % n)
        x,y,d=2,2,1
        while d==1:
            x=g(x,c,n)
            y=g(g(y,c,n),c,n)
            d=m.gcd(abs(x-y),n)
        if d!=n:
            #print(d)
            return d
    return -1

def g(x,c,n):
    return ((x**2)+c) % n

def modular_sqrt(a, p):
    #tonelli shanks algorithm for solving the congruence x^2 = a mod p returns x, -x mod p
    if legendre_symbol(a, p) != 1:
        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return 0
    elif p % 4 == 3:
        x = pow(a, (p + 1) // 4, p)
        return (x, p - x)
    
    s = p - 1
    e = 0
    while s % 2 == 0:
        s //= 2
        e += 1
    
    n = 2
    while legendre_symbol(n, p) != -1:
        n += 1
        
    x = pow(a, (s + 1) // 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e
    
    while True:
        t = b
        m = 0
        for m in range(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return (x, p - x)

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m

def legendre_symbol(a, p):
    #legendre symbol, returns 1 if x^2 = a mod p has a solution, -1 otherwise
    ls = pow(a, (p-1) // 2, p)
    return -1 if ls == p - 1 else ls
    
def isqrt(n):
    #integer square root, finds biggest x such that x^2 <= n
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x


