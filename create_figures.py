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



##Figure 1

n_min,n_max = 0,5
ns = np.arange(n_min,n_max+1)
canalizing_depths = np.arange(max(ns)+1)
proportion_canalizing_depths_all = np.zeros((len(ns),max(ns)+1))

for i,n in enumerate(ns):
    number_f_with_canalizing_depths = [number_k_canalizing_depth(n,k) for k in range(1,n+1)]
    total_fs = 2**(2**n)
    number_f_not_canalizing = total_fs-sum(number_f_with_canalizing_depths)
    proportion_canalizing_depths_all[i,:n+1] = np.append(number_f_not_canalizing,number_f_with_canalizing_depths)/total_fs

fig, ax = plt.subplots(figsize=(5.5,2.5))
for i,canalizing_depth in enumerate(canalizing_depths):
    ax.bar(ns,proportion_canalizing_depths_all[:,i],bottom=np.sum(proportion_canalizing_depths_all[:,:i],1),label=str(canalizing_depth))
ax.legend(frameon=False,loc='center',bbox_to_anchor=[0.5,1.13],ncol=8,title='exact canalizing depth of Boolean function')
ax.set_xticks(ns)
ax.set_xlabel('Number of inputs (n)')
ax.set_ylabel('Proportion of functions')
ax.set_ylim([0,1])
ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'depth_vs_variables_n{n_min}to{n_max}.pdf',bbox_inches='tight')





##Figure 2a

n_max = 5
ns = np.arange(n_max+1)
canalizing_depths = np.arange(max(ns)+1)
proportion_canalizing_depths = np.zeros((len(ns),max(ns)+1))
C_n_k = np.zeros((len(ns),max(ns)+1),dtype=int)
N_n_m_k = np.zeros((len(ns),max(ns)+1,max(ns)+1),dtype=int)

for n in range(n_max+1):
    C_n_k[n,1:n+1] = [number_k_canalizing_depth(n,k) for k in range(1,n+1)]
    total_fs = 2**(2**n)
    C_n_k[n,0] = total_fs-sum(C_n_k[n,1:n+1])
    for m in range(n+1): #number essential variables
        for k in range(m+1): #depth
            if m == k and k == 0:
                N_n_m_k[n,m,k] = 2
            elif k <= m and m < n:
                N_n_m_k[n,m,k] = nchoosek(n,m) * N_n_m_k[m,m,k]
            elif m == n:
                N_n_m_k[n,m,k] = C_n_k[n,k] - sum([N_n_m_k[n,i,k] for i in range(n)])

    proportion_canalizing_depths[n] = N_n_m_k[n,n]/sum(N_n_m_k[n,n])
                
fig, ax = plt.subplots(figsize=(5.5,2.5))
for i,canalizing_depth in enumerate(canalizing_depths):
    ax.bar(ns,proportion_canalizing_depths[:,i],bottom=np.sum(proportion_canalizing_depths[:,:i],1),label=str(canalizing_depth))
ax.legend(frameon=False,loc='center',bbox_to_anchor=[0.5,1.13],ncol=8,title='exact canalizing depth of non-degenerate Boolean function')
ax.set_xticks(ns)
ax.set_xlabel('Number of essential inputs (n)')
ax.set_ylabel('Proportion of functions')
ax.set_ylim([0,1])
ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'depth_vs_essential_variables_n{n_min}to{n_max}.pdf',bbox_inches='tight')





#Figure 2b

width = 0.55
violinplot_args = {'widths': width, 'showmeans': True, 'showextrema': False}
fig, ax = plt.subplots(figsize=(5.5,2.5))

base_gap = 1     # gap between groups
intra_gap = 0.57   # gap within group

max_depth = max(ns)

ax.spines[['right', 'top']].set_visible(False)

positions = []
values = []
colors_used = []
group_centers = []

current_x = 0.0
for i, n in enumerate(ns):
    valid_depths = np.append(np.arange(n-1), n)
    n_viols = len(valid_depths)

    # positions centered on each group's midpoint
    offsets = np.linspace(
        -(n_viols - 1) * intra_gap / 2,
        (n_viols - 1) * intra_gap / 2,
        n_viols
    )
    group_positions = current_x + offsets
    positions.extend(group_positions)
    group_centers.append(current_x)

    for depth in valid_depths:
        values.append(proportion_canalizing_depths[n, depth])
        colors_used.append('C'+str(depth))

    # advance x-position based on total group width
    group_width = (n_viols - 1) * intra_gap
    current_x += group_width / 2 + base_gap + width + intra_gap

    # axis labels
    ax.set_xticks(group_centers)

# plot violins one by one with colors
for vpos, val, c in zip(positions, values, colors_used):
    bar = ax.bar([vpos], val, color=c,width=width)
ax.set_yscale('log')
ax.set_ylabel('Proportion of functions [log]')

# axis labels
ax.set_xlabel('Number of essential inputs (n)')
ax.set_xticks(group_centers)
ax.set_xticklabels(ns)
ax.set_ylim([5e-7,1])
ax.set_yticks([1.e-06, 1.e-04, 1.e-02, 1])
ax.set_yticklabels(['$\\mathdefault{10^{-6}}$','$\\mathdefault{10^{-4}}$','$\\mathdefault{10^{-2}}$','1'])
plt.savefig('proportion_log_nondegenerated.pdf',bbox_inches='tight')






##Figure 3

res = []
for n in ns:
    if n==0:
        continue
    prop_canalizing_all = 1-proportion_canalizing_depths_all[n][0]
    prop_canalizing_nondeg = 1-proportion_canalizing_depths[n][0]
    prop_ncf_all = proportion_canalizing_depths_all[n][n]
    prop_ncf_nondeg = proportion_canalizing_depths[n][n]
    res.append([n,
                prop_canalizing_all,
                prop_canalizing_nondeg,
                (prop_canalizing_all/prop_canalizing_nondeg-1)*100,
                np.log2(prop_canalizing_nondeg/prop_canalizing_all),
                prop_ncf_all,
                prop_ncf_nondeg,
                (prop_ncf_all/prop_ncf_nondeg-1)*100,
                np.log2(prop_ncf_nondeg/prop_ncf_all)
                ])
res = np.array(res)
# np.round(res,4)

# f,ax = plt.subplots()
# ax.plot(res[:,0],res[:,7],'o')
# ax.plot(res[:,0],res[:,3],'x')
# ax.spines[['right', 'top']].set_visible(False)
# ax.set_xlabel('Number of (essential) inputs (n)')
# ax.set_ylabel('Log2-fold change')

# [x1,x2] = ax.get_xlim()
# ax.set_xticks(ns)
# ax.plot([x1,x2],[0,0],'k--')
# ax.set_xlim([x1,x2])






# Compute log2-fold change Δlog2 = log2(tilde_P / P)
log2_fold_canal = np.log2(res[:,1] / res[:,2])
log2_fold_NCF = np.log2(res[:,5] / res[:,6])

f, ax = plt.subplots(figsize=(4.5, 3.0))

ax.plot(res[:,0], log2_fold_canal, 'o-', color='C6', label='canalizing')#label=r'$\log_2(\tilde P / P)$')
ax.plot(res[:,0], log2_fold_NCF, 'x--', color='C9', label='NCF')#label=r'$\log_2(\tilde P / P)$')

# Baseline line at 0 (no bias)
ax.axhline(0, color='k', linestyle='--', linewidth=1)

# Aesthetics
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('Number of (essential) inputs $n$')
ax.set_ylabel(r'$\log_2$-fold change, ' + r'$\log_2(\tilde P_* / P_*)$')
ax.set_xticks(res[:,0])
ax.set_xlim(res[:,0].min() - 0.2, res[:,0].max() + 0.2)

ax.legend(frameon=False)
plt.tight_layout()
plt.savefig('proportion_comp.pdf',bbox_inches='tight')
