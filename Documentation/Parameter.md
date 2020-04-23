## Semiconductor Parameter
| Parameter | Unit | Explanation |
| :------------- |:-------------:| :-----|
| SC_Interface |  | This variable defines the position of interfaces between different layers of semiconductors. The interface is defined in terms of Grid/ Node-Points, i.e. 0<=SC_Interface<=NodeCount_x-1. The first semiconductor interface has to be set to 0.|
| Epsi_s | |          Dielectric constant of the semiconductor |
| Chi |  eV |Electron affinity of the semiconductor |
| m_cb | | Density of states effective mass of electrons in conduction band as a fraction of the electron mass |
| m_vb | | Density of states effective mass of holes in valence band as a fraction of the electron mass |
| E_g | eV | Bandgap of the semiconductor |
| Band_Gap_Type | |Bandgap TYP of the semiconductor (accepted values: Direct/InDirect) |
| N_A | 1/cm^3 | Acceptor concentration |
| N_D | 1/cm^3 | Donator concentration |
| E_D | eV | Energy of donor level (relative to valence band edge) |
| E_A | eV | Energy of acceptor level (relative to valence band edge) |
| mu_n | cm^2/(Vs) | Mobility of electrons |
| mu_p | cm^2/(Vs) | Mobility of holes |
| E_SS | eV | This is the peak position of the gaussian distribution of the surface states (relative to the valence band edge). |
| E_CNL | eV | The charge neutrality level of the surface state distribution. (relative to the valence band edge) |
| FWHM | eV | Full-width-half-maximum value of the gaussian distribution of the surface states |
| surface_charge_density | 1/cm^2 | Surface charge density. This is the area enclosed by the gaussian distribution. This value is in the range between 10^12 (nearly unpinned) and 10^14 (fully pinned). Set this value to 0, if no surface states are required in the calculation. |
| P_x |  1/nm^2 | Defines the polarisation field of the semiconductor in x Direction. |
| P_y | 1/nm^2 | Defines the polarisation field of the semiconductor in y Direction. |
| P_z | 1/nm^2 | Defines the polarisation field of the semiconductor in z Direction. |
| Interface_charge | 1/cm^2 | Interface_charge |
| Surface_charge | 1/cm^2 | Surface_charge |
| alpha | 1/cm | Absorption coefficient |
| tau | s | Minority charrier lifetime |
| light_on | | Switch for turning on the laser (accepted values: true/false) |
| Inv_V | | Allow free holes in the valence band (accepted values: true/false) |
| Inv_C | | Allow free electrons in the conduction band (accepted values: true/false) |
|P_opt | W | Optical power of the laser (set to 0, if no additional excess carriers should be included) |
| E_ph | eV | Energy of a photon (emitted by the laser) |
| A_opt | nm^2 | Illuminated area |

## Simulation Parameter:
| Parameter | Unit | Explanation |
| :------------- |:-------------:| :-----|
| T | K | Temperature (Note T=0 K ist not supported at the moment. Will be fixed soon.)|
| phi_m | eV | Work function of the tip | 
| Tip_Radius | nm | Radius of a hyperbolic tip |
| Tip_Apex_Angle | ° | Opening angle of the tip 0<Angle<90. Smaller value means sharper tip. |
| d | nm | Tip-sample separation |
| NodeCount_x | | Number of grid points in x direction. The minimal Number of grid points is depending on the physical length and the minimal step width. Condition: NodeCount_x> 3*log2(x_length/2/dx). |
| NodeCount_y | | Number of grid points in y direction. The minimal Number of grid points is depending on the physical length and the minimal step width. Condition: NodeCount_y> 3*log2(y_length/2/dy). |
| NodeCount_z | | Number of grid points in z direction (direction normal to the semiconductors surface). The minimal Number of grid points is depending on the physical length and the minimal step width. Condition NodeCount_z> 3*log2(z_length/2/dz). |
| x_length | nm | Physical length of the system in x direction. |
| y_length | nm | Physical length of the system in y direction. |
| z_length | nm | Physical length of the system in z direction. |
| dx_min | nm | Minimal step width of the computation grid in x direction. At the borders of the computation grid the step width will be increased to match the total grid width. |
| dy_min | nm | Minimal step width of the computation grid in y direction. At the borders of the computation grid the step width will be increased to match the total grid width. |
| dz_min | nm |Minimal step width of the computation grid in z direction. At the borders of the computation grid the step width will be increased to match the total grid width. |
| dPhi | | Abort criteria of the numerical solver: If the relative change of the solution in one iteration step is smaller than dPhi, the final result is obtained. (e.g. dPhi = 1e-7) |
| refinement_steps | | Defines the number of grid-node refinements. [1..n] Each refinement doubles the number of grid points in each direction. Be carefull: Each refinement step needs eight times (2^3) the memory of the previous step! |
| stepwidth_z | eV | Stepwidth in the outer integral in Eq. 2 in [2] - Accuracy: 0.04 eV (low quality, may be to large, E_F is near a band edge), 0.01 eV (normal quality). |
| z_stop | nm | Upper limit for z-integration in transmission coefficients. Standard value: 200 nm, large values enhance accuracy of the tunnel currents but take longer. |
| omega | |Relaxation parameter for under-relaxation of the Newton iteration method. For semiconductors, values around 0.2 are suggested by Selberherr[1]. Decrease this value, if the solution diverges. |
| scale_tunnelarea | | Scale the tunnel area (accepted values: true/false) |

##  Output Parameter:
| Parameter | Unit | Explanation |
| :------------- |:-------------:| :-----|
| V_List | V | Voltage list to be calculated seperated by a comma.|
| I(V):output | | Output of the tunneling current I(V) |
| I(V):filename | |Filename for output file of I(V) curve. File is stored in program directory, no subdirectory is specified |
| BB_central(z):output | | band bending as a function of the distance from the semiconductors surface for each voltage step defined above |
| BB_central(z):filename | | Filenames for the band bending output files. Each voltage step produces a new ascii file containing V, E_V(z),E_C(z), E_F, E_F+V. %.2f in the filename will be replaced by the voltage (for more information about the filename see format command in pascal or c language) |
| BB(V):output | | band bending as a function of voltage |
| BB(V):filename | | Filename for the band bending output file. File contains V, E_V(V), E_C(V), E_F + V' |
| EFQ(V):output | | Quasi-Fermi levels and quasi-effective masses as a function of voltage. |
| EFQ(V):filename | | Filename for the Quasi-Fermi level output file |
| Phi3D:output | | output of the whole three dimensional phi-array in binary form (for plotting with gnuplot) and in ASCII-Format. |
| Phi3D:filename | | file extension will be added automatically. (.bin=binary, .txt=ascii). The first parameter (%.2f) is the voltage given as floating point value, The second parameter (%d) is the x-index of the phi(x,y,z) array and given as integer value.|
| Rho3D:output | | output of the whole three dimensional electron and hole concentrations in binary form |
| Rho3D:n_filename | | |
| Rho3D:p_filename | | |
| Rho_central(z):output | | |
| Rho_central(z):filename | | |
| Interim_solution:Intervall | | Saves an iterim solution every k-th iteration step to investigate the evolution of the potential and carrier densities. (0 = no output of iterim solution) |
| Interim_solution:x | |Index (integer value) of the plane, that should be stored (from 0 to NodeCount_x-1). Usually: NodeCount_x/2 |
| Interim_solution:Path | | Folder name for storing the iterim solution of a spezific voltage. The first parameter (%.2f) is the voltage given as floating point value |
| Interim_solution:n_filename | | Filename for electron concentration array file. First parameter (%d) is the iteration step (interger value) |
| Interim_solution:p_filename | | Filename for hole concentration array file. First parameter (%d) is the iteration step (interger value) |
| Interim_solution:phi_filename | | Filename for phi array file. First parameter (%d) is the iteration step (interger value) |
| Interim_solution:ND_filename | | Filename for concentration of ionized donor array file. First parameter (%d) is the iteration step (interger value)|
| Interim_solution:NA_filename | | Filename for concentration of ionized acceptor array file. First parameter (%d) is the iteration step (interger value) |
| Interim_solution:phi_quantum_n_filename | | Filename for quantum potential of electrons array file. First parameter (%d) is the iteration step (interger value) |
| Interim_solution:phi_quantum_p_filename | | Filename for quantum potential of holes array file. First parameter (%d) is the iteration step (interger value) |