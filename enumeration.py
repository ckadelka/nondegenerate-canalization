#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:56:16 2025

@author: ckadelka
"""

import numpy as np
import matplotlib.pyplot as plt


def fak(n):
    if n==0 or n==1:
        return 1
    else:
        return n*fak(n-1)

def ndbinom(n,ks):
    if len(ks)==1:
        return nchoosek(n,ks[0])
    res= 1
    for k in ks:
        res*=fak(k)
    return fak(n)/res

def nbinom(n,k):
    return fak(n)/(fak(k)*fak(n-k))

def Bstar(n):
    return 2**(2**n) - 2*((-1)**n -n) + sum([(-1)**k*nbinom(n,k)*2**(k+1)*2**(2**(n-k)) for k in range(1,n+1)])


def nchoosek(population, sample):
    "Returns `population` choose `sample`."
    s = max(sample, population - sample)
    assert s <= population
    assert population > -1
    if s == population:
        return 1
    numerator = 1
    denominator = 1
    for i in range(s+1, population + 1):
        numerator *= i
        denominator *= (i - s)
    return numerator/denominator

def sterling_times_fak_r(n,r):
    return sum([((-1)**i)*nchoosek(r,i)*(r-i)**n for i in range(r+1)])
    
def sterling_difference(n,r):
    return r**n+sum([(-1)**i*nchoosek(r-1,i-1)*(r-i)**(n-1)*(r**2*1./i-r+n) for i in range(1,r+1)])

def number_ncfs(n,p):
    if n==1 and p==2:
        return (2,0,2)
    if p==2:
        N1=0
    else:
        N1=2**(n-1)*p*(p-2)*(p-1)**n*sum([(p-1)**(r-1)*n*sterling_times_fak_r(n-1,r-1) for r in range(2,n+1)])
    N2=2**n*p*(p-1)**n*sum([(p-1)**r*sterling_difference(n,r) for r in range(1,n)])
    return N1+N2,N1,N2

def partitions(n, I=1):
    yield (n,)
    for i in range(I, n//2 + 1):
        for p in partitions(n-i, i):
            yield (i,) + p

def number_k_canalizing_depth(n,k):
    return nchoosek(n,k)*(number_ncfs(k,2)[0]+Bstar(n-k)*2**(k+1)*sum([ndbinom(k,ks) for ks in partitions(k)])    )

def get_C_n_k(n_max):
    ns = np.arange(n_max+1)
    C_n_k = np.zeros((len(ns),max(ns)+1),dtype=int)    
    for n in range(n_max+1):
        C_n_k[n,1:n+1] = [number_k_canalizing_depth(n,k) for k in range(1,n+1)]
        total_fs = 2**(2**n)
        C_n_k[n,0] = total_fs-sum(C_n_k[n,1:n+1])
    return C_n_k

def get_N_n_m_k(n_max):
    C_n_k = get_C_n_k(n_max)
    ns = np.arange(n_max+1)
    N_n_m_k = np.zeros((len(ns),max(ns)+1,max(ns)+1),dtype=int)
    
    for n in range(n_max+1):
        for m in range(n+1): #number essential variables
            for k in range(m+1): #depth
                if m == k and k == 0:
                    N_n_m_k[n,m,k] = 2
                elif k <= m and m < n:
                    N_n_m_k[n,m,k] = nchoosek(n,m) * N_n_m_k[m,m,k]
                elif m == n:
                    N_n_m_k[n,m,k] = C_n_k[n,k] - sum([N_n_m_k[n,i,k] for i in range(n)])
    return N_n_m_k

def get_proportion_canalizing_depths(n_max,ALLOW_DEGENERATE_FUNCTIONS=False):
    ns = np.arange(n_max+1)
    proportion_canalizing_depths = np.zeros((len(ns),max(ns)+1))
    if ALLOW_DEGENERATE_FUNCTIONS:
        for i,n in enumerate(ns):
            number_f_with_canalizing_depths = [number_k_canalizing_depth(n,k) for k in range(1,n+1)]
            total_fs = 2**(2**n)
            number_f_not_canalizing = total_fs-sum(number_f_with_canalizing_depths)
            proportion_canalizing_depths[n,:n+1] = np.append(number_f_not_canalizing,number_f_with_canalizing_depths)/total_fs
    else:
        N_n_m_k = get_N_n_m_k(n_max)
        for n in range(n_max+1):
            proportion_canalizing_depths[n] = N_n_m_k[n,n]/sum(N_n_m_k[n,n])
    return proportion_canalizing_depths

def get_proportion_canalizing(n_max,ALLOW_DEGENERATE_FUNCTIONS=False):
    ns = np.arange(n_max+1)
    proportion_canalizing_depths = get_proportion_canalizing_depths(n_max,
                                                                    ALLOW_DEGENERATE_FUNCTIONS=ALLOW_DEGENERATE_FUNCTIONS)
    proportion_canalizing = []
    for n in ns:
        if n==0:
            proportion_canalizing.append(np.nan)
            continue
        proportion_canalizing.append(1-proportion_canalizing_depths[n][0])
    return np.array(proportion_canalizing)

def get_proportion_NCF(n_max,ALLOW_DEGENERATE_FUNCTIONS=False):
    ns = np.arange(n_max+1)
    proportion_canalizing_depths = get_proportion_canalizing_depths(n_max,
                                                                    ALLOW_DEGENERATE_FUNCTIONS=ALLOW_DEGENERATE_FUNCTIONS)
    proportion_NCF = []
    for n in ns:
        if n==0:
            proportion_NCF.append(np.nan)
            continue
        proportion_NCF.append(proportion_canalizing_depths[n][n])
    return np.array(proportion_NCF)
