#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:56:16 2025

@author: ckadelka
"""

import numpy as np
import matplotlib.pyplot as plt

import enumeration
              

fig_folder = 'figs/'
n_max = 5
ns = np.arange(n_max+1)

#Computations
proportion_canalizing_depths_all = enumeration.get_proportion_canalizing_depths(n_max,ALLOW_DEGENERATE_FUNCTIONS=True)
proportion_canalizing_depths_nondegenerated = enumeration.get_proportion_canalizing_depths(n_max,ALLOW_DEGENERATE_FUNCTIONS=False)

proportion_canalizing_all = enumeration.get_proportion_canalizing(n_max,ALLOW_DEGENERATE_FUNCTIONS=True)
proportion_canalizing_nondegenerated = enumeration.get_proportion_canalizing(n_max,ALLOW_DEGENERATE_FUNCTIONS=False)

proportion_NCF_all = enumeration.get_proportion_NCF(n_max,ALLOW_DEGENERATE_FUNCTIONS=True)
proportion_NCF_nondegenerated = enumeration.get_proportion_NCF(n_max,ALLOW_DEGENERATE_FUNCTIONS=False)

#log2-fold change Δlog2 = log2(proportion_*_all / proportion_*_nondegenerated)
log2_fold_canalizing = np.log2(proportion_canalizing_all / proportion_canalizing_nondegenerated)
log2_fold_NCF = np.log2(proportion_NCF_all / proportion_NCF_nondegenerated)


##Figure 1
fig, ax = plt.subplots(figsize=(5.5,2.5))
for i,canalizing_depth in enumerate(ns):
    ax.bar(ns,proportion_canalizing_depths_all[:,i],bottom=np.sum(proportion_canalizing_depths_all[:,:i],1),label=str(canalizing_depth))
ax.legend(frameon=False,loc='center',bbox_to_anchor=[0.5,1.13],ncol=8,title='exact canalizing depth of Boolean function')
ax.set_xticks(ns)
ax.set_xlabel('Number of inputs (n)')
ax.set_ylabel('Proportion of functions')
ax.set_ylim([0,1])
ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{fig_folder}depth_vs_variables_n0to{n_max}.pdf',bbox_inches='tight')


##Figure 2a
fig, ax = plt.subplots(figsize=(5.5,2.5))
for i,canalizing_depth in enumerate(ns):
    ax.bar(ns,
           proportion_canalizing_depths_nondegenerated[:,i],
           bottom=np.sum(proportion_canalizing_depths_nondegenerated[:,:i],1),
           label=str(canalizing_depth))
ax.legend(frameon=False,loc='center',bbox_to_anchor=[0.5,1.13],ncol=8,title='exact canalizing depth of non-degenerate Boolean function')
ax.set_xticks(ns)
ax.set_xlabel('Number of essential inputs (n)')
ax.set_ylabel('Proportion of functions')
ax.set_ylim([0,1])
ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{fig_folder}depth_vs_essential_variables_n0to{n_max}.pdf',bbox_inches='tight')


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
        values.append(proportion_canalizing_depths_nondegenerated[n, depth])
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
plt.savefig(f'{fig_folder}proportion_log_nondegenerated_n0to{n_max}.pdf',bbox_inches='tight')


##Figure 3



f, ax = plt.subplots(figsize=(4.5, 3.0))

ax.plot(ns, log2_fold_canalizing, 'o-', color='C6', label='canalizing')#label=r'$\log_2(\tilde P / P)$')
ax.plot(ns, log2_fold_NCF, 'x--', color='C9', label='NCF')#label=r'$\log_2(\tilde P / P)$')

# Baseline line at 0 (no bias)
ax.axhline(0, color='k', linestyle='--', linewidth=1)

# Aesthetics
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('Number of (essential) inputs $n$')
ax.set_ylabel(r'$\log_2$-fold change, ' + r'$\log_2(\tilde P_* / P_*)$')
ax.set_xticks(ns)
ax.set_xlim(ns[np.isnan(log2_fold_canalizing)==False].min() - 0.2, 
            ns[np.isnan(log2_fold_canalizing)==False].max() + 0.2)

ax.legend(frameon=False)
plt.tight_layout()
plt.savefig(f'{fig_folder}proportion_comp.pdf',bbox_inches='tight')
