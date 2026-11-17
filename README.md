# APPROXIMATION OF THE LEVY-DRIVEN STOCHASTIC HEAT EQUATION ON THE SPHERE
Simulation and visualization of Lévy random field and stochastic heat equation on the unit sphere <i>S^2</i> and convergence rate of the spectral and Euler-Maruyama approximation.
This repository represents a snapshot of the code used for the creation of the manuscript

<i>APPROXIMATION OF THE LEVY-DRIVEN STOCHASTIC HEAT EQUATION ON THE SPHERE</i>

It allows simulating sample paths of Lèvy random field and stochastic heat equation in time using Spectral methods in space and Euler-Maruyama in time.
Further, it contains functions for computing the convergence rate of the methods either direct or with Monte-Carlo simulations.
The script are in MATLAB and they generate plots (or videos) with the mentioned results. 

# Content of the repository

levy.m containts the function to select the type of noise. Different noise possibilities are explained directly in the file.
It is needed to have this function in the same folder as the other files.

Regarding the other files, detailed explanations are provided in the comments within each .m file. However, for convenience, here is a brief overview of what each script does:

<i>LevyProcessOnSphere_Movie.m </i><br>
Generate movies and images of samples of the Levy process on the sphere for a specific time, grid, noise and initial condition.

<i>Samples_stochastic_heat_equation_and_driving_noise.m </i><br>
Generate movies and images of samples of the Levy process and the corresponding sample of the stochastic heat equation on the sphere for a specific time, grid, noise and initial condition.
The movies and images are created using an exponential transformation on the sphere.

<i>SpectralApprox_StrongError_MeanError_SecondMomentError.m</i><br>
Compute and generate plots for the convergence rate of the spectral approximation for the strong error, the mean error and the second moment error.

<i>SpectralApprox_WeakError.m</i><br>
Compute and generate plots for the weak convergence rate using Monte Carlo sampling of the backward Euler-Maruyama method for the simulation and test functions of the form |u(x)|^p where |.| represents the L^2 norm on the sphere.

<i>TemporalApprox_Euler_Direkt_Strong_Mean_SecondMoment.m</i><br>
Compute and generate plots for the convergence rate of the temporal approximation using direct simulation method for the backward Euler-Maruyama scheme for the strong error, the mean error and the second moment error.

# Requirements
This code was run under MATLAB R2024b.

# Authors and acknowledgment
This code was created by Annika Lang, Andrea Papini and Verena Schwarz.
This research was funded in parts by the Austrian Science Fund (FWF) [10.55776/DOC78], by the European Union (ERC, StochMan, 101088589), by the Swedish Research Council (VR) through grant no. 2020-04170, by the Wallenberg AI, Autonomous Systems and Software Program (WASP) funded by the Knut and Alice Wallenberg Foundation, and by the Chalmers AI Research Center (CHAIR). For open access purposes, the authors have applied a CC BY public copyright license to any author-accepted manuscript version arising from this submission. 

# License
This snapshot is published under the GPL v3.0 license.

# Disclaimer
Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them. For open access purposes, the authors have applied a CC BY public copyright license to any author-accepted manuscript version arising from this submission. 
