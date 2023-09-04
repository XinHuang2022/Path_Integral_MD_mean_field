# Implementation of path integral molecular dynamics approximation

## Purpose

This repository contains the code for implementation of the path integral method with Feynman-Kac formulation, to approximate the quantum canonical partition function $Z$ and mean-field $h$ for a fermionic system consisting of indistinguishable electrons.
Two cases with different forms of potential are tested, namely:
* `Case V1`: Confining harmonic oscillator potential $V(\mathbf{x})=\frac{1}{2}|\mathbf{x}|^2$, and,
* `Case V2`: Coulomb interactions between electrons under an external harmoic trap, $V(\mathbf{x})= \frac{1}{2}|\mathbf{x}|^2+\frac{1}{2}\sum_{k}\sum_{l\ne k}\frac{\lambda}{|x_k-x_l|}$.

## Getting Started

The main function can be accessed through the file named `TesterPIMC_BB.m`, with the following specification of the input arguments:
*  `dim0`: the dimensionality of the space.
*  `Ne0`: the number of electrons.
*  `beta0`: the inverse temperature $1 / ( k_B * T )$.
*  `dt0`: the time step for the integration in the exponent of the Feymann-Kac representation.
*  `Mxmax0`: an integer for the exponent of the sample size $N$ for the independent estimators of partition function $Z$ and mean-field $h$, i.e., $N$ equals to $2$ to the power of $\mathrm{Mxmax0}$.
*  `dMx0`: when plotting the figures to test the convergence of estimators with increasing sample size, the distance between two data points in the $x$-axis.
*  `Nestim0`: the number of independent replicas of partition function $Z$ and mean-field $h$ estimators, with each replica generates and stores 2^Mxmax0 samples to approximate the Monte Carlo integral.
*  `Ncore0`: the number of parallel workers for evaluating all the $2$^$\mathrm{Nestim0}$ estimators.

A simple example with implementation at the inverse temperature $\beta=1$ with $n=6$ particles in dimension $d=3$ utilising $4$ CPU-cores can be achieved with the following input values:
| dim0 | Ne0 | beta0 | dt0 | Mxmax0 | dMx0 | Nestim0 | Ncore0 |
|----------|----------|----------|----------|----------|----------|----------|----------|
| 3 | 6 | 1 | 0.025 | 12 | 0.5 | 8 | 4 |


With $\mathrm{Mxmax0} = 12$ and $\mathrm{Nestim0}=8$, a total sample set of size $2^{12+8}=2^{20}$ of estimators will be generated, and then the function `stat_and_plot.m` shall collect all the generated samples and compute the sample mean value, and also give the statistical confidence intervals for the mean-field estimator $h$.

## Usage Reference




 
