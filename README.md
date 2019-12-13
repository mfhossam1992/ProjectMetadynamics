# ProjectMetadynamics

This is a tookit written in C++, which performs hybrid Metadynamics - Molecular Dynamics Simulation, using biased potential. 
This allows one exploring the Energy Landscape of the system at hand, as a function of a set of chosen collective variables, within the framework of Metadynamics. 

The potential Energy Surface is biased by adding a set of gaussians, centered, in the collective variable space, at different values of the collective variable, chosen/and/recorded each tau time-steps. This history-dependent gaussians allow one to reconstruct the free energy surface as a function of the chosen biasing collective variable, by inverting this history-dependent gaussian summed surface.

This code was developed by me, Hossam Farag, along with couple of other colleagues, Navi Ning and Harry Feldman, at UIUC, as part of a class project for MSE-485.

