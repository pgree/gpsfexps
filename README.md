# Generalized prolate spheroidal functions (GPSFs) subroutines

This repository includes Fortran codes for computing with generalized prolate 
spheroidal functions (GPSFs), which are the extension of prolate spheroidal 
wave functions to greater than one dimension. These codes are an implementation of 
the algorithms of https://arxiv.org/abs/1811.02733 and include subroutines 
for evaluating GPSFs, quadrature of bandlimited functions, and approximation 
via GPSF expansion. 

Codes for computing with GPSFs are included in `gpsfexps.f` and were
written by Philip Greengard and Kirill Serkh. 
All `*.f` files other than `gpsfexps.f` and `xk_exps.f` were written by Vladimir 
Rokhlin who generously provided them for use in this repository.

As of 10/9/2022, these codes have not been seriously optimized. 
There are several potential speedups that could be made including use of 
fast numerical linear algebra software. 

These codes were originally written and tested in a Linux environment using
the "f90" compiler and have also been tested on several Macs using the 
gfortran compiler. 

## Numerical tests of equations in the paper

In `gpsfexps.f` we include tests of a few equations in the paper. We reference
formulas using the original arxiv version of the paper https://arxiv.org/abs/1811.02733.

- Evaluating GPSFs: the subroutine `plot_phi_88` plots the radial component of GPSFs
  defined in (87) and (88). Plotting is done using gnuplot.
- Eigenfunctions of linear operator $L_{N, c}$: subroutine `test_lnc_93` checks equation (93)
- Eigenfunctions of integral operator $M_{N, c}$: subroutine `test_mnc_90` checks equation (90)
- Evaluation of eigenvalues of GPSF integral operator: subroutine `test_gammas_4` checks the algorithms
  of Sections 4.1 and 4.2
- Second-order ODE satisfied by GPSFs: subroutine `test_ode_236` checks equation (236)
  
