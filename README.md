# Implementation of path integral molecular dynamics approximation

## Purpose

This repository contains the code for implementation of the path integral method with Feynman-Kac formulation, to approximate the quantum canonical partition function $Z$ and mean-field h for a fermionic system consisting of indistinguishable electrons.
Two cases with different forms of potential are tested, namely:
*'Case V1': Confining harmonic oscillator potential $V(\mathbf{x})=\frac{1}{2}|\mathbf{x}|^2$, and,
*'Case V2': Coulomb interactions between electrons under an external harmoic trap, $V(\mathbf{x})= \frac{1}{2}|\mathbf{x}|^2+\frac{1}{2}\sum_{k}\sum_{l\ne k}\frac{\lambda}{|x_k-x_l|}$
under the harmonic oscillator potential, for multiple number of fermionic or distinguishable particles.


 
