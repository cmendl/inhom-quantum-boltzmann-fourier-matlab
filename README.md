Simulation of a spatially inhomogeneous matrix-valued quantum Boltzmann equation using Fourier transformation
=============================================================================================================

Matlab implementation, test files and visualization routines for the simulations in "Numerical scheme for a spatially inhomogeneous matrix-valued quantum Boltzmann equation" (see References). The main simulation files are `simHomBoltzmann.m` (spatially homogeneous) and `simInhomBoltzmann.m` (spatially inhomogeneous) in the *code* subfolder.

The *visualization* subfolder contains routines for visualizing the simulation output files. Call `loadHomMatlab('../output/simdata_hom_N32_example.mat')` and `loadInhomMatlab('../output/simdata_inhom_N32_Dirichlet_example.mat')` for the spatially homogeneous and inhomogeneous simulation examples, respectively.

Note that an alternative C implementation with the same functionality is also available online.


License
-------
Copyright (c) 2014, Christian B. Mendl  
All rights reserved.  
http://christian.mendl.net

This program is free software; you can redistribute it and/or
modify it under the terms of the Simplified BSD License
http://www.opensource.org/licenses/bsd-license.php


References
----------
1. Jianfeng Lu, Christian B. Mendl  
   Numerical scheme for a spatially inhomogeneous matrix-valued quantum Boltzmann equation  
   Journal of Computational Physics 291, 303-316 (2015), [arXiv:1408.1782](http://arxiv.org/abs/1408.1782)
2. Martin L.R. Fürst, Christian B. Mendl, Herbert Spohn  
   Matrix-valued Boltzmann equation for the Hubbard chain  
   Physical Review E 86, 031122 (2012), [arXiv:1207.6926](http://arxiv.org/abs/1207.6926)
