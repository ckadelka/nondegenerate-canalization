# Non-degenerate Canalizing Boolean Functions

This repository contains Python code to reproduce all figures and numerical results from the paper  
**"On the number of non-degenerate canalizing Boolean functions"** by Claus Kadelka.

## Overview
The scripts compute and visualize:
- Recursive enumeration of Boolean functions by number of essential variables and canalizing depth (`N(n, m, k)`).
- Prevalence of canalizing and nested canalizing functions among all and non-degenerate functions.
- Log2-fold change quantifying the bias introduced when degeneracy is ignored.

## Contents
- `enumeration.py` – recursive computation of `N(n, m, k)` values.  
- `figures.py` – plotting scripts for all figures in the manuscript.  
- `data/` – precomputed results tables (optional).  
- `notebooks/` – interactive Jupyter notebooks reproducing results.

## Requirements
Python ≥ 3.9 with `numpy`, `matplotlib`, and `pandas`.

## Citation
If you use this code, please cite:
> Claus Kadelka (2025). *On the number of non-degenerate canalizing Boolean functions.* Physica D (submitted).  
> DOI: [link to arXiv preprint]

## License
Released under the MIT License. See `LICENSE` for details.
