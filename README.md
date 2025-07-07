# APPROXIMATION OF THE LEVY-DRIVEN STOCHASTIC HEAT EQUATION ON THE SPHERE
Simulation and visualization of Lévy random field and stochastic heat equation on the unit sphere <i>S^2</i> and convergence rate of the spectral and Euler-Maruyama approximation.
This repository represents a snapshot of the code used for the creation of the manuscript

<i>APPROXIMATION OF THE LEVY-DRIVEN STOCHASTIC HEAT EQUATION ON THE SPHERE</i>

It allows simulating sample paths of Lèvy random field and stochastic heat equation in time using Spectral methods in space and Euler-Maruyama in time.
Further, it contains functions for computing the convergence rate of the methods either direct or with Monte-Carlo simulations.
The script are in MATLAB and they generate plots (or videos) with the mentioned results. 

# Content of the repository

Levy.m containts the function to select the type of noise. We have different (12) noise possibilities explained directly in the file.
It is needed to have this function in the same folder as the other files.
Regarding the other files, everything is written in comments inside of each ".m"-file, but as a small rundown here it is what each code does:

<i>Samples_stochastic_heat_equation_and_driving_noise.m </i><br>
Generate movies and snapshots of the Levy Random Fields and the Stochastic Heat equation for a specific time, grid, noise and initial condition.
The movie are created using an exponential transformation on the sphere, i.e. x in S^2 is moved to x*exp(u(x)/|u(x)|), where u is the solution to the equation.

<i>Strong_Mean_SecondMoment_Convergence_Rate_Spectral_SHE_Levy.m</i><br>
Compute and generate plot for the convergence rate of the spectral scheme for a reference value of the truncation series.
In particular using the4 parameter ExpLP, VarLP is it possible to compute the Strong rate (1,1), the Mean rate (1,0) and Second Moment (0,1).

<i>Spectral_Weak_Error_Rate_SHE_Independent_Levy.m</i><br>
For a given test function of the form |u(x)|^p where |.| represents the L^2 norm on the sphere, we compute the weak rate of convergence for the spectral method.
The computation is done using Monte-Carlo sampling (M=20, changable in the file), with a reference solution at fixed time grid and with different truncation parameters.
PorAselector is a variable that can swwitch between observing the decay for different test function (i.e. different p, select the variable to 0) and observing the decay for a fixed test fanction (with selected p) respect to the regolarity of the field (i.e. different alpha, select the variable to 1).

<i>StrongSpectralError_ForwardEM_BackwardEM_Direct.m</i><br>
Compute the strong error (mean and second moment) for the Euler-Maruyama scheme, using direct computation for the mean and variance and with both the forward and backward scheme.
The scheme are sequenced, so for each run 3 plots will be shown: Strong Spectral Error, Forward EM Strong Error, Backward EM Strong Error.

<i>StrongSpectralError_MonteCarloForwardEM_MonteCarloBackwardEM.m</i><br>
Compute the strong error (mean and second moment) for the Euler-Maruyama scheme, using Monte-Carlo sampling for the mean and variance and with both the forward and backward scheme.
The scheme are sequenced, so for each run 3 plots will be shown: Strong Spectral Error, Forward EM Strong Error, Backward EM Strong Error.

# Requirements
This code was run under MATLAB R2024b.

# Authors and acknowledgment
This code was created by Annika Lang, Andrea Papini and Verena Schwarz.
This research was funded in parts by the Austrian Science Fund (FWF) [10.55776/DOC78], by the European Union (ERC, StochMan, 101088589), by the Swedish Research Council (VR) through grant no. 2020-04170, by the Wallenberg AI, Autonomous Systems and Software Program (WASP) funded by the Knut and Alice Wallenberg Foundation, and by the Chalmers AI Research Center (CHAIR). For open access purposes, the authors have applied a CC BY public copyright license to any author-accepted manuscript version arising from this submission. 

# License
This snapshot is published under the GPL v3.0 license.

# Disclaimer
Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them. For open access purposes, the authors have applied a CC BY public copyright license to any author-accepted manuscript version arising from this submission. 
