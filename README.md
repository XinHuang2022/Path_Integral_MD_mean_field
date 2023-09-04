# Implementation of path integral molecular dynamics

## Purpose

This repository contains the **MATLAB** code for implementation of the path integral method with Feynman-Kac formulation, to approximate the quantum canonical partition function $Z$ and mean-field $h$ for a fermionic system consisting of indistinguishable electrons.
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
*  `Nrep0`: the number of independent replicas of partition function $Z$ and mean-field $h$ estimators, with each replica generates and stores $2^\mathrm{Mxmax0}$ samples to approximate the Monte Carlo integral.
*  `Ncore0`: the number of parallel workers for evaluating all the $2^{\mathrm{Nrep0}}\times 2^\mathrm{Mxmax0}$ estimators.

A simple example with implementation at the inverse temperature $\beta=1$ with $n=6$ particles in dimension $d=3$ can be achieved with the following input values, utilising $4$ CPU-cores.
| dim0 | Ne0 | beta0 | dt0 | Mxmax0 | dMx0 | Nrep0 | Ncore0 |
|----------|----------|----------|----------|----------|----------|----------|----------|
| 3 | 6 | 1 | 0.025 | 12 | 0.5 | 8 | 4 |


With $\mathrm{Mxmax0} = 12$ and $\mathrm{Nrep0}=8$, a total sample set of size $2^{12+8}=2^{20}$ of estimators will be generated, and then the function `stat_and_plot.m` shall collect all the generated samples and compute the sample mean value, and also give the statistical confidence intervals for the mean-field estimator $h$.

## Usage Reference

### Regarding the choice of sample size
Since the Monte Carlo method is applied to approximate a high-dimensional integral, the sample size parameter $\mathrm{Mxmax0}$ should be sufficiently large.

Specifically in the code, it is assumed that $\mathrm{Mxmax0}\geq 6$, so that at least $2^6=64$ independent estimators are generated for each replica. Otherwise if $\mathrm{Mxmax0}< 6$ the code may return an error, and the statistics part will not be implemented. Moreover, as the inverse temperature parameter $\beta$ increases, the sample size should accordingly increase fast, due to the exponentially increasing sample variance.

### Regarding the choice of time step size $\Delta t$
In this implementation for simplicity, the number of time slices $M=\beta/\Delta t$ on the region $[0,\beta]$ is assumed to be an integer. 

For example, given $\beta=1.5$, selecting $\Delta t=0.03$ leads to $M=50$ time slices, which produces a good approximation. 

However for $\beta=1$ with $\Delta t=0.03$ we have $M=33.3$ which is non-integer and can lead to further errors. It is recommended to always choose $\Delta t$ such that $M$ is an integer, hence for the latter case, a recommended time step value is $\Delta t = 0.025$ or $\Delta t = 0.0125$, leading to $M=40$ or $M=80$ time slices respectively.

### Regarding the parallelization with the number of replicas
In order to achieve an optimized efficiency, it is recommended to have $2^{\mathrm{Nrep0}}$ larger than the number of CPU-cores available for this computation, so that the workload between all the processors are balanced.

Besides, if this code is run with **MATLAB** version `R2023a`, the largest number of parallel workers is recommended to not exceed 100, i.e., $\mathrm{Ncore0}\leq 100$.  This problem is hopefully fixed with the updated version of **MATLAB** `R2023b`.




 
