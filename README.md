# GUI for the Simulation *P_SpaceChargeLight*

**Introduction**

P_SpaceChargeLight is a solver for both, Poisson- and continuity equations for holes and electrons within a semiconductor[1]. Like the program SEMITIP of R. M. Feenstra, P_SpaceChargeLight is designed especially for the solution of problems involving a metalic probe tip in proximity to a semiconductor surface. Nevertheless, P_SpaceChargeLight was developed "from scratch" and has no source code in common with SEMITIP. Due to the solution of the continuity equations, P_SpaceChargeLight is capable of calculating the distribution of additional (excess) carriers, that are influenced by the electric field of the tip. The carrier concentrations within the semiconductor are the effective-mass (or parabolic) approximation. The numerical Poisson- and continuity equation solver is based upon a finite-difference scheme in three dimensions. Rotational symmetries are not included due to reasons of universality (e. g. non-symmetrically tips) The tunnel currents are computed in Bardeens approximation[2], which was further developed and adapted by Harrison[3], Bono and Good[4], as well as Feenstra[5]. The potential and carrier concentrations along the central axis through the tip-apex is used to derive the tunnel current in a one dimensional approximation.  Excess carriers are included in the tunnel currents by the derivation of quasi Fermi levels.

Further information about the physical and numerical concepts is available in the technical manual.

This software is freely provided for the use by any researcher. Users are encouraged to report bugs and suggestions for improvements and to specify the software ("P_SpaceChargeLight", with reference to this website and to the publication [1]) in their publications.

**Running the program**

P_SpaceChargeLight can be downloaded as source-code and precompiled executable for Microsoft Windows 32/64 Bit, Linux, and MacOS. The source-code is written in the Pascal language. It can be compiled with the freely available Free-Pascal-Compiler (FPC) and also with the commercially available Embarcadero Delphi/ RAD Studio environment.

The program P_SpaceCharge needs the input file "parameters.txt" to be present in the execution directory. This ASCII-File contains all relevant information and physical parameters of the semiconductor (bulk and surface), the tip, and the simulation grid. A precise description of all parameters in the "parameters.txt" file is given here.

**Physical Background and Theory**

The theoretical models can be found in the references. Additionally, a precise description of the models, finite-difference iteration schemes, and interface- and boundary conditions is available as [PDF](https://iffgit.fz-juelich.de/mschned/p_spacechargelight/-/blob/master/Documentation/The_Physics_Of_P_SpaceChargeLight.pdf). Special attention should be paid to the book of Prof. S. Selberherr[6] for the derivation of finite-difference and finite-element equations for semiconductors. Especially the implementation of the boundary conditions is often done in an erroneous way, leading to increased numiercal errors at the boundaries. Due to very small tip-radii and thus small tip-induced quantum dots, the classical Drift-and-Diffusion model that is used to evaluate the carrier current densities may overestimate the carrier densities and thus one has to take into account a quantum correction: P_SpaceChargeLight derives additional repulsive potentials for holes and electrons to account for the quantum compressibility[7]. The additional potentials are then used as a perturbation in the continuity equations.

**Output files**

P_SpaceChargeLight is able to generate a variety of different output files like potential distributions, electron- and hole distributions, band bending along the central axis, Quasi-Fermi levels, and tunnel currents. It depends on the settings within the parameter.txt file, which kind of data is stored. (see the description of the parameter file for further information.) Each output of data can be switched on or off, and individual filenames can be assigned. A detailed description of the individual files can be found here.

**References**


[1] M. Schnedler, V. Portz, P. H. Weidlich, R. E. Dunin-Borkowski, and Ph. Ebert, *Quantitative description of photoexcited scanning tunneling spectroscopy and its application to the GaAs(110) surface*, Phys. Rev. B **91**, 235305

[2] J. Bardeen, *Tunnelling from a many-particle point of view*, Phys. Rev. Lett. **6**,  57-59

[3] W. A. Harrison, *Tunneling from an Independent-Particle Point of View*, Phys. Rev. **123**, 85

[4] J. Bono and R. H. Good, *Theoretical discussion of the scanning tunneling microscope applied to a semiconductor surface*, Surf. Sci. **175**, 415-420

[5] R. M. Feenstra and J. A. Stroscio, *Tunneling spectroscopy of the GaAs(110) surface*, J. Vac. Sci. Technol. B **5**, 923-929

[6] S. Selberherr, *Analysis and Simulation of Semiconductor Devices*, Springer Vienna-New York, ISBN 978-3-7091-8754-8, 1984

[7] N. Rowsey, R. P. Muller, R. P. Young, *3D TCAD Modeling of Candidate Structures for the Silicon Qubit*, CSRI Summer Proceedings 2009 (Computer Science Research Institute at Sandia National Laboratories, Albuquerque, New Mexico, USA)
