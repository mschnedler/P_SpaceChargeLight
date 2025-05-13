unit U_Constants;

interface
const
       ProgramVersion = '10/18/2024 (V27, Extensions: Light, Quantum, Polar, FDMSolve, OhmicContact)';
       {$IFDEF FPC}
        TargetOS = {$I %FPCTARGETOS%};
       {$ELSE}
        {$IFDEF WIN32}
          TargetOS = 'Windows 32 bit';
        {$ENDIF}
        {$IFDEF WIN64}
          TargetOS = 'Windows 64 bit';
        {$ENDIF}
        {$IFDEF MACOS}
          TargetOS = 'MacOS X 32/64 bit';
        {$ENDIF}
       {$ENDIF}

       Integration_Accuracy=1e-12;
       kb     = 1.38064852e-23;      // Boltzmann constant in J/K
       k      = 8.61733035e-05;        // Bolzmann constant in eV/K
       e      = 1.60217662e-19;       // Electron charge in C
       h      = 6.62607004e-34;    // Planck constant in J*s  //vorher: 6.6260755e-34
       h_bar  = h /(2*Pi);        // reduced Planck constant in J*s
       me     = 9.10938e-31;      // Restmass of electron in kg
       Epsi_0 = 8.85419e-12;      // Permittivity of vacuum in As/(V*m)
       Epsi_vac = Epsi_0 * 1e-9;  // Vacuum permittivity in [C/(V*nm)]
       eVToJ  = 1.6022e-19;       // Factor between eV and Joule in Joule/eV
       newton_steps_max = 100;    // Steps used in Newton algorithm to find zero. (e.g. in Find_EF)
       Precision        = 0.0000001; // First non-significant digit in the iteration result
       h_eV   = h / eVToJ;        // Planck constant in eV*s

       FDI_Coefficients: array[0..3] of double =
        //(3.53553e-1,   -4.95009e-3,   1.48386e-4,   -4.42563e-6);
       (3.53553390e-1,-4.95008973e-3,1.48385771e-4,-4.42563012e-6);

       ParameterNames : array[0..75] of string = ('T','d','phi_m','epsi_s','chi',
                                                 'm_cb','m_vb','E_g','N_D','N_A',
                                                 'E_D','E_A','Inv_V','Inv_C','P_x','P_y','P_z',
                                                 'stepwidth_z',
                                                 'I(V):output','V_list',
                                                 'I(V):filename','BB_central(Z):output',
                                                 'BB_central(Z):filename','BB(V):output',
                                                 'BB(V):filename','Phi3D:output',
                                                 'Phi3D:filename','NodeCount_x',
                                                 'NodeCount_y', 'NodeCount_z',
                                                  'x_length','y_length','z_length',
                                                  'dx_min','dy_min','dz_min',
                                                  'dPhi','omega',
                                                  'Tip_Radius','Tip_Apex_Angle',
                                                  'refinement_steps','E_SS','E_CNL','FWHM',
                                                  'surface_charge_density','z_stop','mu_n',
                                                  'mu_p','P_opt','E_ph','A_opt','alpha',
                                                  'tau_n','tau_p','light_on',
                                                  'Rho3D:output','Rho3D:n_filename',
                                                  'Rho3D:p_filename','Rho_central(z):output',
                                                  'Rho_central(z):filename',
                                                  'Interim_solution:Intervall',
                                                  'Interim_solution:x',
                                                  'Interim_solution:Path',
                                                  'Interim_solution:n_filename',
                                                  'Interim_solution:p_filename',
                                                  'Interim_solution:phi_filename',
                                                  'SC_Interface','Band_Gap_Type',
                                                  'Interim_solution:ND_filename',
                                                  'Interim_solution:NA_filename',
                                                  'Interim_solution:phi_quantum_n_filename',
                                                  'Interim_solution:phi_quantum_p_filename',
                                                  'EFQ(V):output','EFQ(V):filename',
												                          'Interface_charge','Surface_charge');
       ErrorMessages: array[0..4] of string = ('Number of array elements differs from previous declaration.',
                                               'Newton iteration diverged.','Index out of CB-DOS list.',
                                               'Index out of VB-DOS list.','Not a valid value for Band_Gap_Type: Direct/Indirect expected'
                                              );
implementation
end.


