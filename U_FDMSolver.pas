unit U_FDMSolver;
// ********************************************************************
// * Unit U_FDMSolver.pas                                             *
// *                                                                  *
// * Self-consistent evaluation of Poisson- and Continuity equations  *
// * in semiconductors and at semiconductor/vacuum/metal interfaces   *
// *                                                                  *
// * Author:                                                          *
// *                                                                  *
// * Michael Schnedler                                                *
// * Forschungszentrum Juelich GmbH                                   *
// * Wilhelm-Johnen-Strasse                                           *
// * 52428 Juelich                                                    *
// *                                                                  *
// * Please cite: 10.1103/PhysRevB.91.235305                          *
// *              10.1103/PhysRevB.93.195444                          *
// *                                                                  *
// * Published under GNU General Public License v3.0                  *
// ********************************************************************
interface

uses
       U_Constants,
       U_Common,
       U_CarrierConcentration,
       Math,
       System.DateUtils,
       SysUtils;

const
       FDMSolver_Version = '4 (06/28/2024)';

type
       TGridPointArray=array of double;
       TProgressNotification = procedure(k,m:integer;CheckPoint_Phi,relChange,T:double);stdcall;

var
       checkpoint_x,checkpoint_y,checkpoint_z:integer;  // Definition of position x,y,z of Phi value which is displayed for every iteration step
       SpaceChargeArray: TSpaceChargeArray;             // Array that stores all 3D information of our system
       NodeCount_x, NodeCount_y, NodeCount_z: integer;  // Number of Points in x,y,z direction of the space charge grid
       x_length,y_length, z_length:           double;   // Width of the calculation area in m

       dPhi_Grenz:                            double;   // Limit of accuracy in the phi iteration (relative limit)
       px_array,py_array,pz_array:            TGridPointArray; // Physical positions of node points in x,y, and z-direction
       omega:                                 double;    // Relaxation parameter for Newton SOR iteration method
                                                         // This factor should be set to values around ]0.2[ (according to Selberherr)
       Nscale :                               double;    // Scaling factor for carrier concentrations (for numerical stability of computation) [1/nm^3]
       PhiScale:                              double;    // Scaling factor the electrostatic potential (for numerical stability of computation) [1/V]

       ProgressNotification: TProgressNotification;

       function Phi(x,y,z:integer):double;
       function Get_N(x,y,z:integer):double;
       function Get_P(x,y,z:integer):double;
       procedure Iterate;
       procedure Initialize_Semiconductors;
       procedure SetNodeCount(NodeCountX,NodeCountY,NodeCountZ:integer);
       procedure AssignSC(x,y,z:integer;SC:integer);
       procedure AssignVacuum(x,y,z:integer);
       procedure AssignMetal(x,y,z:integer;Phi:double);
       procedure AssignOhmic(x,y,z:integer;SC:integer; Phi:double);

       // Bibliography
       //  [1] Seiwatz and Green, J. Appl. Phys. 29 (7), 1958
       //  [2] Selberherr, "Analysis and Simulation of Semiconductor Devices",
       //      ISBN 9783709187548, Springer (Wien, New York), 1984
       //  [3] R. M. Feenstra, Semitip V6 source code,
       //      http://www.andrew.cmu.edu/user/feenstra/semitip_v6/
       //  [4] E. Fred Schubert, "Light-Emitting Diodes",
       //      Second Edition, Cambridge University Press, 2006
       //      ISBN 1139455222, 9781139455220
       //  [5] D. Z. Garbuzov, "Reradiation effects, lifetimes and probabilities of
       //      band-to-band transitions in direct A3B5 compounds of GaAs type"
       //      Journal of Luminescence 27, 109-112 (1982)
       //      ISSN 0022-2313
       //      http://dx.doi.org/10.1016/0022-2313(82)90033-3
       //  [6] W. B. Joyce and R. W. Dixon, Appl. Phys.Lett. 31, pp. 354 (1977)

       // Definition of functions and procedures


       // Usage:   1. Call initialize_semiconductors
       //          2. Call SetNodeCount for initialization of px,py,pz arrays and SpaceChargeArray
       //          3. Manually set px_array, py_array, pz_array values
       //          4. Manually set SpaceChargeArray values

implementation



     // Derivation of Diffusion coefficient for electrons under impurity
     // scattering, carrier-carrier scattering, and velocity saturation.
     function Dn(x,y,z:integer):double;
      var beta, vs, mu,mu_c,Ex,Ey,Ez:double;
          SC:integer;
          CI:double;
          Cn_ref,Sn:double;
      begin
          //Comment out the following two lines to
          //activate the spatial dependence of the diffusion coefficient:
          result:=1;
          exit;

          SC:=SpaceChargeArray[x,y,z].SCIndex;
          mu:=Semiconductors[SC].mu_n;

          //ionized impurity scattering Eq. 4.1-23 from Ref [2] p. 87
          CI:=abs(ND_ionized(GetE_Fqn2(SC,Phi(x,y,z), Get_N(x,y,z)),Phi(x,y,z),SC))
             +abs(NA_ionized(GetE_Fqp2(SC,Phi(x,y,z), Get_P(x,y,z)),Phi(x,y,z),SC));
          Cn_ref:=3e16*1e-21; //nm^-3
          Sn:=350;
          mu:=mu/sqrt(1+(CI/(Cn_ref+CI/Sn)));

          // carrier-carrier scattering Eq. 4.1-31 from Ref [2] p. 90
          mu_c:=1.428e20*1e-7/(sqrt(Get_N(x,y,z)*Get_P(x,y,z))*
                ln(1+4.54e11*1e-14*power(Get_N(x,y,z)*Get_P(x,y,z),1/3)));
          mu_c:=mu_c*1e-18;
          mu:=1/(1/mu+1/mu_c);

          // velocity saturation Eq. 4.1-48 from Ref [2] p. 95:
          beta:=2;
          vs:= 2.4e7/(1+0.8*exp(T/600))*1e7; // nm/s
          if x=0 then
            Ex:= (Phi(x+1,y,z)-Phi(x,y,z))/(px_array[x+1]-px_array[x])
          else if x=NodeCount_x-1 then
            Ex:= (Phi(x,y,z)-Phi(x-1,y,z))/(px_array[x]-px_array[x-1])
          else
            Ex:= (Phi(x+1,y,z)-Phi(x-1,y,z))/(px_array[x+1]-px_array[x-1]);

          if y=0 then
            Ey:= (Phi(x,y+1,z)-Phi(x,y,z))/(py_array[y+1]-py_array[y])
          else if y=NodeCount_y-1 then
            Ey:= (Phi(x,y,z)-Phi(x,y-1,z))/(py_array[y]-py_array[y-1])
          else
            Ey:= (Phi(x,y+1,z)-Phi(x,y-1,z))/(py_array[y+1]-py_array[y-1]);

          if z=0 then
            Ez:= (Phi(x,y,z+1)-Phi(x,y,z))/(pz_array[z+1]-pz_array[z])
          else if z=NodeCount_z-1 then
            Ez:= (Phi(x,y,z)-Phi(x,y,z-1))/(pz_array[z]-pz_array[z-1])
          else
            Ez:= (Phi(x,y,z+1)-Phi(x,y,z-1))/(pz_array[z+1]-pz_array[z-1]);

          result:=kb*T/e*1e18*mu/power(1+power(1e18*mu*sqrt(Ex*Ex+Ey*Ey+Ez*Ez)/vs,beta),1/beta);
          result:=result/Semiconductors[SC].Dn;
      end;


     // Derivation of Diffusion coefficient for holes under impurity
     // scattering, carrier-carrier scattering, and velocity saturation.
     function Dp(x,y,z:integer):double;
     var beta, vs, mu,Ex,Ey,Ez:double;
          SC:integer;
          CI,Cp_ref,Sp:Double;
          mu_c:double;
     begin
          //Comment out the following two lines to
          //activate the spatial dependence of the diffusion coefficient:
          result:=1;
          exit;

          SC:=SpaceChargeArray[x,y,z].SCIndex;
          mu:=Semiconductors[SC].mu_p;

          //ionized impurity scattering Eq. 4.1-23 from Ref [2] p. 87
          CI:=abs(ND_ionized(GetE_Fqn2(SC,Phi(x,y,z), Get_N(x,y,z)),Phi(x,y,z),SC))
             +abs(NA_ionized(GetE_Fqp2(SC,Phi(x,y,z), Get_P(x,y,z)),Phi(x,y,z),SC));

          Cp_ref:=4e16*1e-21; //nm^-3
          Sp:=81;
          mu:=mu/sqrt(1+(CI/(Cp_ref+CI/Sp)));


          // carrier-carrier scattering Eq. 4.1-31 from Ref [2] p. 90
          mu_c:=1.428e20*1e-7/(sqrt(Get_N(x,y,z)*Get_P(x,y,z))*
                ln(1+4.54e11*1e-14*power(Get_N(x,y,z)*Get_P(x,y,z),1/3)));
          mu_c:=mu_c*1e-18;

          mu:=1/(1/mu+1/mu_c);

          // velocity saturation Eq. 4.1-48 from Ref [2] p. 95:
          beta:=1;
          vs:= 2.4e7/(1+0.8*exp(T/600))*1e7; // nm/s
          if x=0 then
            Ex:= (Phi(x+1,y,z)-Phi(x,y,z))/(px_array[x+1]-px_array[x])
          else if x=NodeCount_x-1 then
            Ex:= (Phi(x,y,z)-Phi(x-1,y,z))/(px_array[x]-px_array[x-1])
          else
            Ex:= (Phi(x+1,y,z)-Phi(x-1,y,z))/(px_array[x+1]-px_array[x-1]);

          if y=0 then
            Ey:= (Phi(x,y+1,z)-Phi(x,y,z))/(py_array[y+1]-py_array[y])
          else if y=NodeCount_y-1 then
            Ey:= (Phi(x,y,z)-Phi(x,y-1,z))/(py_array[y]-py_array[y-1])
          else
            Ey:= (Phi(x,y+1,z)-Phi(x,y-1,z))/(py_array[y+1]-py_array[y-1]);

          if z=0 then
            Ez:= (Phi(x,y,z+1)-Phi(x,y,z))/(pz_array[z+1]-pz_array[z])
          else if z=NodeCount_z-1 then
            Ez:= (Phi(x,y,z)-Phi(x,y,z-1))/(pz_array[z]-pz_array[z-1])
          else
            Ez:= (Phi(x,y,z+1)-Phi(x,y,z-1))/(pz_array[z+1]-pz_array[z-1]);

          result:=kb*T/e*1e18*mu/power(1+power(1e18*mu*sqrt(Ex*Ex+Ey*Ey+Ez*Ez)/vs,beta),1/beta);
          result:=result/Semiconductors[SC].Dp;
     end;


      procedure Update_Diffusioncoefficients(x,y,z:integer);
      begin
          SpaceChargeArray[x,y,z].DiffN:=Dn(x,y,z);
          SpaceChargeArray[x,y,z].DiffP:=Dp(x,y,z);
      end;

      procedure AssignSC(x,y,z:integer;SC:integer);
      begin
         SpaceChargeArray[x,y,z].Material:=mSC;
         SpaceChargeArray[x,y,z].SCIndex:=SC;
         SpaceChargeArray[x,y,z].Phi:=-(Semiconductors[SC].E_offset);
         SpaceChargeArray[x,y,z].dPhi:=0;
         SpaceChargeArray[x,y,z].n:=Semiconductors[SC].n0;
         SpaceChargeArray[x,y,z].dn:=0;
         SpaceChargeArray[x,y,z].p:=Semiconductors[SC].p0;
         SpaceChargeArray[x,y,z].dp:=0;
         SpaceChargeArray[x,y,z].Phi_Gamma_n:=0;
         SpaceChargeArray[x,y,z].Phi_Gamma_p:=0;
         SpaceChargeArray[x,y,z].Polarisation.x:= Semiconductors[SC].Polarisation.x;
         SpaceChargeArray[x,y,z].Polarisation.y:= Semiconductors[SC].Polarisation.y;
         SpaceChargeArray[x,y,z].Polarisation.z:= Semiconductors[SC].Polarisation.z;
         SpaceChargeArray[x,y,z].OhmicContact:=False;
      end;

     procedure AssignOhmic(x,y,z:integer;SC:integer; Phi:double);
     var C,nisq:double;
     begin
          SpaceChargeArray[x,y,z].OhmicContact:=True;
          // Boundary condition of Ohmic Contact:
          // Carrier concentrations are assumed to be in thermal equilibrium.
          // See Eq. 5.1-12 / 5.1-13 in [2], which is, in our case,
          // equal to:
          SpaceChargeArray[x,y,z].Phi:=-(Semiconductors[SC].E_offset+Phi);
          C:=ND_ionized(GetE_Fqn2(SC,SpaceChargeArray[x,y,z].Phi+Semiconductors[SC].E_offset,SpaceChargeArray[x,y,z].n),SpaceChargeArray[x,y,z].Phi+Semiconductors[SC].E_offset,SC)
            -NA_ionized(GetE_Fqp2(SC,SpaceChargeArray[x,y,z].Phi+Semiconductors[SC].E_offset,SpaceChargeArray[x,y,z].p),SpaceChargeArray[x,y,z].Phi+Semiconductors[SC].E_offset,SC);
          nisq:=Semiconductors[SC].n0*Semiconductors[SC].p0;
          SpaceChargeArray[x,y,z].n:=0.5*(sqrt(C*C+4*nisq)+C);
          SpaceChargeArray[x,y,z].p:=0.5*(sqrt(C*C+4*nisq)-C);

     end;

      procedure AssignVacuum(x,y,z:integer);
      begin
         SpaceChargeArray[x,y,z].Material:=mVac;
         SpaceChargeArray[x,y,z].Phi:=0;
         SpaceChargeArray[x,y,z].n:=0;
         SpaceChargeArray[x,y,z].p:=0;
         SpaceChargeArray[x,y,z].dPhi:=0;
         SpaceChargeArray[x,y,z].dn:=0;
         SpaceChargeArray[x,y,z].dp:=0;
         SpaceChargeArray[x,y,z].Phi_Gamma_n:=0;
         SpaceChargeArray[x,y,z].Phi_Gamma_p:=0;
         SpaceChargeArray[x,y,z].Polarisation.x:= 0;
         SpaceChargeArray[x,y,z].Polarisation.y:= 0;
         SpaceChargeArray[x,y,z].Polarisation.z:= 0;
      end;

      procedure AssignMetal(x,y,z:integer;Phi:double);
      begin
         SpaceChargeArray[x,y,z].Material:=mMetal;
         SpaceChargeArray[x,y,z].Phi:=Phi;
         SpaceChargeArray[x,y,z].n:=0;
         SpaceChargeArray[x,y,z].p:=0;
         SpaceChargeArray[x,y,z].dPhi:=0;
         SpaceChargeArray[x,y,z].dn:=0;
         SpaceChargeArray[x,y,z].dp:=0;
         SpaceChargeArray[x,y,z].Phi_Gamma_n:=0;
         SpaceChargeArray[x,y,z].Phi_Gamma_p:=0;
         SpaceChargeArray[x,y,z].Polarisation.x:= 0;
         SpaceChargeArray[x,y,z].Polarisation.y:= 0;
         SpaceChargeArray[x,y,z].Polarisation.z:= 0;
      end;

      procedure SetNodeCount(NodeCountX,NodeCountY,NodeCountZ:integer);
      begin
        NodeCount_x:=NodeCountX;
        NodeCount_y:=NodeCountY;
        NodeCount_z:=NodeCountZ;
        setlength(SpaceChargeArray,NodeCount_x,NodeCount_y,NodeCount_z);
        setlength(px_array,NodeCount_x);
        setlength(py_array,NodeCount_y);
        setlength(pz_array,NodeCount_z);
      end;


      procedure SetScaling(n_max,p_max,phi_max:double);
      begin
         C_kT := 1/(k*T);
         Nscale:=1/max(n_max,p_max);
         PhiScale:= e/(kb*T);
      end;

      // Calculates the Fermi energies, carrier concentrations and futher
      // properties of all semiconductors.
      // Returns the highest carrier concentrations n_max and p_max of all
      // semiconductors. These values are used for scaling the solution.
      // For details about the scaling see main procedure of this program and [2]
      procedure Initialize_Semiconductors;
      var i           :integer;
          c           :double;
          n_max,p_max,phi_max:double;
      begin
          kT:=k*T;
          Init_FjList(-5,5,10000);
          output('\nInitializing '+inttostr(Semiconductor_Count)+' semiconductors:\n');
          c:=2*Pi*me*kb*T/(h*h)*1e-18; //[1/nm^2]
          n_max:=0;
          p_max:=0;
          phi_max:=0;
          for i := 0 to Semiconductor_Count-1 do
            with Semiconductors[i] do
               begin
                  C_ND:=N_D*1e-27;
                  C_NA:=N_A*1e-27;
                  N_C:=2*power(Semiconductors[i].m_cb*c,3/2);
                  N_V:=2*power(Semiconductors[i].m_vb*c,3/2);
                  if light_on then
                      G_photo:=alpha*P_opt/(E_ph*A_opt) // [1/(nm^3*s)]
                  else
                      G_photo:=0;
                  Semiconductors[i].E_F:=Find_EF(i);
                  n0:=n(E_F,0,inv_c,i);
                  p0:=p(E_F,0,inv_v,i);
                  n_max:=max(n_max,C_ND);
                  p_max:=max(p_max,C_NA);
                  Dn:= mu_n*kb*T/e*1e18;
                  Dp:= mu_p*kb*T/e*1e18;
                  // C_Rate is the Bimolecular Recombination Coefficient
                  // For a direct-band gap semiconductor, the following equation
                  // is used, as proposed by [4] and [5].
                  if (BandGapType=bgtDirect) and (T>0) then
                       C_Rate:= 3e-10*power(300/T,3/2)*power(E_G/1.5,2) * 1e21 //in  [nm^3/s]
                  else
                       C_Rate:=0;
                  //Auger recombination coefficients:
                  c_auger_n:=2.8e-31 * 1e42;  // in [nm^6/s]
                  c_auger_p:=0.99e-31 * 1e42; // in [nm^6/s]
                  //Shockley Reed Hall recombination coefficients:
                  nref_n:=3e17*1e-21; // [nm^-3] reference doping concentration
                  nref_p:=3e17*1e-21; // [nm^-3] reference doping concentration
                  tau_p:=tau_p0/(1+(N_D+N_A)*1e-27/nref_p);
                  tau_n:=tau_n0/(1+(N_D+N_A)*1e-27/nref_n);

                  C_Poisson:=e/(Epsi_0*Epsi_s)*1e9;
                  Epsi_semi :=Epsi_0*Epsi_s*1e-9;
                  Init_SurfRhoList(10000,i);
               end;


          for i := 0 to Semiconductor_Count-1 do
            begin
              Semiconductors[i].E_offset:= Semiconductors[i].E_g
                                          -Semiconductors[i].E_f
                                          +Semiconductors[i].Chi
                                          -(Semiconductors[0].E_g
                                            -Semiconductors[0].E_f
                                            +Semiconductors[0].Chi);

              Semiconductors[i].Discontinuity_VB:= (Semiconductors[0].Chi+Semiconductors[0].E_g)

                                                  -(Semiconductors[i].Chi+Semiconductors[i].E_g);

              Semiconductors[i].Discontinuity_CB:= Semiconductors[0].Chi-Semiconductors[i].Chi;

              Phi_max:=max(phi_max,abs(Semiconductors[i].E_offset));

            end;

          SetScaling(n_max,p_max,phi_max);
      end;


      // Returns the divergence of the polarisation at the integer position x,y,z.
      // Unit is 1/nm^3 (devided by e)
      procedure Prepare_Div_P;
      var dPdx,dPdy,dPdz:double;
          x,y,z:integer;
      begin
        for x := 0 to NodeCount_x-1 do
          for y := 0 to NodeCount_y-1 do
            for z := 0 to NodeCount_z-1 do
              if SpaceChargeArray[x,y,z].Material=mSC then
                begin
                  if (x=0) or ((x>0) and (x<NodeCount_x-1) and (SpaceChargeArray[x-1,y,z].SCIndex<>SpaceChargeArray[x,y,z].SCIndex)) then
                    dPdx:=(SpaceChargeArray[x+1,y,z].Polarisation.x-SpaceChargeArray[x,y,z].Polarisation.x)/(px_array[x+1]-px_array[x])
                  else if (x=NodeCount_x-1) or ((x>0) and (x<NodeCount_x-1) and (SpaceChargeArray[x+1,y,z].SCIndex<>SpaceChargeArray[x,y,z].SCIndex)) then
                    dPdx:=(SpaceChargeArray[x,y,z].Polarisation.x-SpaceChargeArray[x-1,y,z].Polarisation.x)/(px_array[x]-px_array[x-1])
                  else
                    dPdx:=(SpaceChargeArray[x+1,y,z].Polarisation.x-SpaceChargeArray[x-1,y,z].Polarisation.x)/(px_array[x+1]-px_array[x-1]);
                  if y=0 then
                    dPdy:=(SpaceChargeArray[x,y+1,z].Polarisation.y-SpaceChargeArray[x,y,z].Polarisation.y)/(py_array[y+1]-py_array[y])
                  else if y=NodeCount_y-1 then
                    dPdy:=(SpaceChargeArray[x,y,z].Polarisation.y-SpaceChargeArray[x,y-1,z].Polarisation.y)/(py_array[y]-py_array[y-1])
                  else
                    dPdy:=(SpaceChargeArray[x,y+1,z].Polarisation.y-SpaceChargeArray[x,y-1,z].Polarisation.y)/(py_array[y+1]-py_array[y-1]);
                  if (z=0) or ((z>0) and (z<NodeCount_z-1) and (SpaceChargeArray[x,y,z-1].Material=mVac)) then
                    dPdz:=(SpaceChargeArray[x,y,z+1].Polarisation.z-SpaceChargeArray[x,y,z].Polarisation.z)/(pz_array[z+1]-pz_array[z])
                  else if (z=NodeCount_z-1) or ((z>0) and (z<NodeCount_z-1) and (SpaceChargeArray[x,y,z+1].Material=mVac)) then
                    dPdz:=(SpaceChargeArray[x,y,z].Polarisation.z-SpaceChargeArray[x,y,z-1].Polarisation.z)/(pz_array[z]-pz_array[z-1])
                  else
                    dPdz:=(SpaceChargeArray[x,y,z+1].Polarisation.z-SpaceChargeArray[x,y,z-1].Polarisation.z)/(pz_array[z+1]-pz_array[z-1]);
                  SpaceChargeArray[x,y,z].divP:=dPdx+dPdy+dPdz;
                end
                else SpaceChargeArray[x,y,z].divP:=0;
      end;


      // Returns the electrostatic potential Phi at the integer position x,y,z.
      // (to be more precise: it is the electrostatic potential of the previous
      //  iteration step: k-1)
      function Phi(x,y,z:integer):double;
      begin
          result:=SpaceChargeArray[x,y,z].Phi;
      end;

      // Returns the electrostatic potential Phi at the integer position x,y,z.
      // (to be more precise: it is the electrostatic potential of the current
      //  iteration step: k)
      function Phi_new(x,y,z:integer):double;
      begin
          result:=SpaceChargeArray[x,y,z].Phi+SpaceChargeArray[x,y,z].dPhi;
      end;

      // Returns the electron concentration n at the integer position x,y,z.
      // (to be more precise: it is the electron concentration of the previous
      //  iteration step: k-1)
      function Get_N(x,y,z:integer):double;
      begin
            result:=SpaceChargeArray[x,y,z].n;
      end;

      // Returns the electron concentration n at the integer position x,y,z.
      // (to be more precise: it is the electron concentration of the current
      //  iteration step: k)
      function Get_N_new(x,y,z:integer):double;
      begin
          result:=SpaceChargeArray[x,y,z].n+SpaceChargeArray[x,y,z].dn;
      end;

      // Returns the hole concentration p at the integer position x,y,z.
      // (to be more precise: it is the hole concentration of the previous
      //  iteration step: k-1)
      function Get_P(x,y,z:integer):double;
      begin
            result:=SpaceChargeArray[x,y,z].p;
      end;

      // Returns the hole concentration p at the integer position x,y,z.
      // (to be more precise: it is the hole concentration of the current
      //  iteration step: k)
      function Get_P_new(x,y,z:integer):double;
      begin
          result:=SpaceChargeArray[x,y,z].p+SpaceChargeArray[x,y,z].dp;
      end;

      function Gamma_n(SC:integer; Phi:double; n:double):double;
      const
          C=1.1283791671; // C=2/sqrt(Pi)
      var E_Fqn:double;
          eta_n:double;
          a,b:double;
      begin
         E_Fqn:=GetE_Fqn(SC,Phi,n);
         eta_n:=(E_Fqn-(Semiconductors[SC].E_g-Phi))/kT;
         if eta_n<-40 then
           begin
             result:=1;
             exit;
           end;
         a:=GetFromFjList(eta_n*kT);
         b:=exp(-eta_n);
         result:=a*b*C;
      end;

      function Gamma_p(SC:integer; Phi:double; p:double):double;
      const
          C=1.1283791671; // C=2/sqrt(Pi)
      var E_Fqp:double;
          eta_p:double;
          a,b:double;
      begin
         E_Fqp:=GetE_Fqp(SC,Phi,p);
         eta_p:=((-Phi)-E_Fqp)/kT;
         if eta_p<-40 then
           begin
             result:=1;
             exit;
           end;
         a:=GetFromFjList(eta_p*kT);
         b:=exp(-eta_p);
         result:=a*b*C;
      end;

      function Phi_Gamma_n(x,y,z:integer; UsePhiNew:boolean = false):double;
      begin
        result:=-Semiconductors[SpaceChargeArray[x,y,z].SCIndex].Discontinuity_CB;
        result:=result+SpaceChargeArray[x,y,z].Phi_Gamma_n;
        result:=result+SpaceChargeArray[x,y,z].Phi;
        if UsePhiNew then
          result:=result+SpaceChargeArray[x,y,z].dPhi;
      end;

      function Phi_Gamma_p(x,y,z:integer; UsePhiNew:boolean = false):double;
      begin
        result:=-Semiconductors[SpaceChargeArray[x,y,z].SCIndex].Discontinuity_VB;
        result:=result+SpaceChargeArray[x,y,z].Phi_Gamma_p;
        result:=result+SpaceChargeArray[x,y,z].Phi;
        if UsePhiNew then
          result:=result+SpaceChargeArray[x,y,z].dPhi;
      end;

      function Update_Phi_Gamma_n(x,y,z:integer):double;
      var SC:integer;
      begin
          SC:=SpaceChargeArray[x,y,z].SCIndex;
          result:=kT*ln(Gamma_n(SpaceChargeArray[x,y,z].SCIndex,
                                SpaceChargeArray[x,y,z].Phi+Semiconductors[SC].E_Offset,
                                SpaceChargeArray[x,y,z].n));
      end;

      function Update_Phi_Gamma_p(x,y,z:integer):double;
      var SC:integer;
      begin
          SC:=SpaceChargeArray[x,y,z].SCIndex;
          result:=-kT*ln(Gamma_p(SpaceChargeArray[x,y,z].SCIndex,
                                 SpaceChargeArray[x,y,z].Phi+Semiconductors[SC].E_Offset,
                                 SpaceChargeArray[x,y,z].p));
      end;



      procedure DetectSurfacesAndInterfaces;
      var x,y,z:integer;
          SurfCount,IntCount:integer;
          SurfX,SurfY,SurfZ,IntX,IntY,IntZ:boolean;
      begin
        for x := 0 to NodeCount_x-1 do
          for y := 0 to NodeCount_y-1 do
            for z := 0 to NodeCount_z-1 do
              begin
                SurfCount:=0;
                IntCount:=0;
                SurfX:= (x>0) and (x<NodeCount_x-1) and (SpaceChargeArray[x,y,z].Material=mSC);
                IntX:=SurfX;
                SurfX:=SurfX  and
                        ((SpaceChargeArray[x-1,y,z].Material=mVac) or
                         (SpaceChargeArray[x+1,y,z].Material=mVac));
                IntX:=IntX and
                        (SpaceChargeArray[x-1,y,z].Material=mSC) and
                        (SpaceChargeArray[x-1,y,z].SCIndex<>SpaceChargeArray[x,y,z].SCIndex);

                SurfY:= (y>0) and (y<NodeCount_y-1) and (SpaceChargeArray[x,y,z].Material=mSC);
                IntY:=SurfY;
                SurfY:=SurfY  and
                        ((SpaceChargeArray[x,y-1,z].Material=mVac) or
                         (SpaceChargeArray[x,y+1,z].Material=mVac));
                IntY:=IntY and
                        (SpaceChargeArray[x,y-1,z].Material=mSC) and
                        (SpaceChargeArray[x,y-1,z].SCIndex<>SpaceChargeArray[x,y,z].SCIndex);

                SurfZ:= (z>0) and (z<NodeCount_z-1) and (SpaceChargeArray[x,y,z].Material=mSC);
                IntZ:=SurfZ;
                SurfZ:=SurfZ  and
                        ((SpaceChargeArray[x,y,z-1].Material=mVac) or
                         (SpaceChargeArray[x,y,z+1].Material=mVac));
                IntZ:=IntZ and
                        (SpaceChargeArray[x,y,z-1].Material=mSC) and
                        (SpaceChargeArray[x,y,z-1].SCIndex<>SpaceChargeArray[x,y,z].SCIndex);

                if SurfX then inc(SurfCount);
                if SurfY then inc(SurfCount);
                if SurfZ then inc(SurfCount);
                if IntX then inc(IntCount);
                if IntY then inc(IntCount);
                if IntZ then inc(IntCount);
                SpaceChargeArray[x,y,z].SurfX:=SurfX;
                SpaceChargeArray[x,y,z].SurfY:=SurfY;
                SpaceChargeArray[x,y,z].SurfZ:=SurfZ;
                SpaceChargeArray[x,y,z].IntX:=IntX;
                SpaceChargeArray[x,y,z].IntY:=IntY;
                SpaceChargeArray[x,y,z].IntZ:=IntZ;
                SpaceChargeArray[x,y,z].SurfCount:=SurfCount;
                SpaceChargeArray[x,y,z].IntCount:=IntCount;
               // if IntCount>0 then
               //    setlength(SpaceChargeArray[x,y,z].Int_n_p,6);
              end;

      end;



      // Helper function for F1 that includes modifications of the Poission
      // equation at surfaces or interfaces of the semiconductor.

      // Three different situations are treated in each spatial direction:

      // 1. Nodepoint is a surface
      // 2. Nodepoint is an interface (e.g. at a hetero-/ homojunction)
      // 3. Nodepoint is a surface and an interface

      // 1. Case 1: "Surface"
      // If the Nodepoint is at a surface only, the discretiation of the
      // Poission equation is modified in analogy to Eq. 6.1-130 in Ref [2].
      // A surface charge Q_int is calculated by assuming a Gaussian distribution
      // of the surface state density, which is integrated from the charge-
      // neutrality level to the Fermi-level (=Q_surf). A change in Polarization
      // is accounted for, too.

      // 2. Case 2: "Interface":
      // Interfaces are treated similar to surfaces, but with two
      // fundamental differences: I. The permittivity to the left and to the
      // right side of the interface are given by the respective semiconductors.
      // II. In Eq. 6.1-129 in Ref [2], a respective charge term is added. This
      // yields that rhos of both semiconductors (left and right to the interface)
      // need to be calculated at the position of the interface (i.e. at x,y,z).
      // While the charge of the "right" semiconductor is stored in n,p which
      // are parameters of the function F1, the charge of the "left" semiconductor
      // is stored separately in Int_n_p[0..6] = n,p,Na,Nd, dNa_dPhi, dNd_dPhi.
      // A change in Polarization is accounted for, too.

      // 3. "Multiple surfaces and interfaces"
      // If there are more than one interface and surface at a
      // single nodepoint, the following approach is used:
      // The boundary condition n*(D1-D2) + Qint = 0 (cf. Eq. 6.1-125 in Ref [2])
      // where n is the normal vector and D1, D2 are the dielectric displacement
      // fields on both sides of the interface and Qint is the interface/
      // surface charge has to be fullfilled for all interfaces/surfaces.
      // The dielectric displacement fields are substituted by the first
      // derivative of the electrostatic potential phi in normal direction n.
      // For the first derivative, only first order one-sided terms are
      // used, while the second order terms, as shown in Eq. 6.1-126 in Ref [2],
      // are ignored.

      procedure SurfaceAndInterfaces(x,y,z:integer; aPhi,n,p, E_Offset,
                                    px,py,pz,pxp1,pxm1,pyp1,pym1,pzp1,pzm1:double;
                                    SC,SurfCount,IntCount:integer;
                                    var PartX,PartY,PartZ, rho: double);
      var SurfX,SurfY,SurfZ,IntX,IntY,IntZ:boolean;
          Cfrac:double;
          Epsi1,Epsi2:double;
          Q_int, Phi_quasi, P_diff:double;
          Q_surf:double;
          IShift:integer;
      begin

        if SurfCount+IntCount>1 then
          begin
            PartX:=0;
            PartY:=0;
            PartZ:=0;
          end;

        SurfX:=SpaceChargeArray[x,y,z].SurfX;
        SurfY:=SpaceChargeArray[x,y,z].SurfY;
        SurfZ:=SpaceChargeArray[x,y,z].SurfZ;
        IntX:=SpaceChargeArray[x,y,z].IntX;
        IntY:=SpaceChargeArray[x,y,z].IntY;
        IntZ:=SpaceChargeArray[x,y,z].IntZ;

        P_diff:=0;
        Q_surf:=0;
        Cfrac:=0;

        // Determination of surface state induced surface charge
        if SurfCount>0 then
        begin
           phi_quasi:=Semiconductors[SC].E_F-GetE_Fqn(SC,aPhi+E_offset,n)-(aPhi+E_offset);
           if not GetRhoSurf(phi_quasi,Q_surf,SC) then output('out of rhosurflist.');
           Q_surf:=(e*Semiconductors[SC].Const_Surface_charge+Q_surf)*1e-18;
        end;

        if SurfX or IntX then
        begin
           // Determination of Permittivity and Polarization change on both
           // sides of the interface/surface
           Epsi1:=Epsi_vac;
           Epsi2:=Epsi_vac;
           if SpaceChargeArray[x+1,y,z].Material=mVac then IShift:=1 else IShift:=0;
           if SpaceChargeArray[x+IShift,y,z].Material =mSC then
               Epsi1:= Semiconductors[SpaceChargeArray[x+IShift,y,z].SCIndex].Epsi_semi;
           if SpaceChargeArray[x-1+IShift,y,z].Material =mSC then
               Epsi2:= Semiconductors[SpaceChargeArray[x-1+IShift,y,z].SCIndex].Epsi_semi;
           P_diff:=-e*(SpaceChargeArray[x+IShift,y,z].Polarisation.x-
                       SpaceChargeArray[x-1+IShift,y,z].Polarisation.x);

           // Assignment of the surface charge density
           if SurfX then Q_int:=Q_surf else Q_int:=0;

           if SurfCount+IntCount=1 then
             begin
               // Case 1 and 2: Surface or Interface
               PartX:= ( Epsi1*(Phi(x+1,y,z)-aPhi)/(pxp1-px)
                         -Epsi2*(aPhi-Phi(x-1,y,z))/(px-pxm1)
                         +Q_int+P_diff
                         )/((Epsi1*(pxp1-px)+Epsi2*(px-pxm1))/2);
               if SurfX then
                 if SpaceChargeArray[x-1+IShift,y,z].Material=mVac then
                    Cfrac:=Epsi1*(pxp1-px)/(Epsi1*(pxp1-px)+Epsi2*(px-pxm1))
                 else
                    Cfrac:=Epsi2*(px-pxm1)/(Epsi1*(pxp1-px)+Epsi2*(px-pxm1));
               if IntX then
                 begin
                    rho:= (pxp1-px)*rho + (px-pxm1)*rho;
                    rho:=rho/(Epsi1*(pxp1-px)+Epsi2*(px-pxm1))*e;
                    rho:=rho/Semiconductors[SC].C_Poisson;
                    Cfrac:=1;
                 end;
             end
           else
             // Case 3: Multiple surfaces and interfaces
             PartX:=   Epsi1*(Phi(x+1,y,z)-aPhi)/(pxp1-px)
                       -Epsi2*(aPhi-Phi(x-1,y,z))/(px-pxm1)
                       +Q_int+P_diff;
        end;

        if SurfY or IntY then
        begin
           // Determination of Permittivity and Polarization change on both
           // sides of the interface/surface
           Epsi1:=Epsi_vac;
           Epsi2:=Epsi_vac;
           if SpaceChargeArray[x,y+1,z].Material=mVac then IShift:=1 else IShift:=0;
           if SpaceChargeArray[x,y+IShift,z].Material =mSC then
               Epsi1:= Semiconductors[SpaceChargeArray[x,y+IShift,z].SCIndex].Epsi_semi;
           if SpaceChargeArray[x,y-1+IShift,z].Material =mSC then
               Epsi2:= Semiconductors[SpaceChargeArray[x,y-1+IShift,z].SCIndex].Epsi_semi;
           P_diff:=-e*(SpaceChargeArray[x,y+IShift,z].Polarisation.y-
                       SpaceChargeArray[x,y-1+IShift,z].Polarisation.y);

           // Assignment of the surface charge density
           if SurfY then Q_int:=Q_surf else Q_int:=0;

           if SurfCount+IntCount=1 then
             begin
               // Case 1 and 2: Surface or Interface
               PartY:= ( Epsi1*(Phi(x,y+1,z)-aPhi)/(pyp1-py)
                         -Epsi2*(aPhi-Phi(x,y-1,z))/(py-pym1)
                         +Q_int+P_diff
                         )/((Epsi1*(pyp1-py)+Epsi2*(py-pym1))/2);
               if SurfY then
                 if SpaceChargeArray[x,y-1+IShift,z].Material=mVac then
                    Cfrac:=Epsi1*(pyp1-py)/(Epsi1*(pyp1-py)+Epsi2*(py-pym1))
                 else
                    Cfrac:=Epsi2*(py-pym1)/(Epsi1*(pyp1-py)+Epsi2*(py-pym1));
               if IntY then
                 begin
                    rho:= (pyp1-py)*rho+(py-pym1)*rho;
                    rho:=rho/(Epsi1*(pyp1-py)+Epsi2*(py-pym1))*e;
                    rho:=rho/Semiconductors[SC].C_Poisson;
                    Cfrac:=1;
                 end;
             end
           else
             // Case 3: Multiple surfaces and interfaces
             PartY:=   Epsi1*(Phi(x,y+1,z)-aPhi)/(pyp1-py)
                       -Epsi2*(aPhi-Phi(x,y-1,z))/(py-pym1)
                       +Q_int+P_diff;
        end;

        if SurfZ or IntZ then
        begin
           // Determination of Permittivity and Polarization change on both
           // sides of the interface/surface
           Epsi1:=Epsi_vac;
           Epsi2:=Epsi_vac;
           if SpaceChargeArray[x,y,z+1].Material=mVac then IShift:=1 else IShift:=0;
           if SpaceChargeArray[x,y,z+IShift].Material =mSC then
               Epsi1:= Semiconductors[SpaceChargeArray[x,y,z+IShift].SCIndex].Epsi_semi;
           if SpaceChargeArray[x,y,z-1+IShift].Material =mSC then
               Epsi2:= Semiconductors[SpaceChargeArray[x,y,z-1+IShift].SCIndex].Epsi_semi;
           P_diff:=-e*(SpaceChargeArray[x,y,z+IShift].Polarisation.z-
                       SpaceChargeArray[x,y,z-1+IShift].Polarisation.z);

           // Assignment of the surface charge density
           if SurfZ then Q_int:=Q_surf else Q_int:=0;

           if SurfCount+IntCount=1 then
             begin
               // Case 1 and 2: Surface or Interface
               PartZ:= ( Epsi1*(Phi(x,y,z+1)-aPhi)/(pzp1-pz)
                         -Epsi2*(aPhi-Phi(x,y,z-1))/(pz-pzm1)
                         +Q_int+P_diff
                         )/((Epsi1*(pzp1-pz)+Epsi2*(pz-pzm1))/2);
               if SurfZ then
                 if SpaceChargeArray[x,y,z-1+IShift].Material=mVac then
                    Cfrac:=Epsi1*(pzp1-pz)/(Epsi1*(pzp1-pz)+Epsi2*(pz-pzm1))
                 else
                    Cfrac:=Epsi2*(pz-pzm1)/(Epsi1*(pzp1-pz)+Epsi2*(pz-pzm1));
               if IntZ then
                 begin
                    rho:= (pzp1-pz)*rho+(pz-pzm1)*rho;
                    rho:=rho/(Epsi1*(pzp1-pz)+Epsi2*(pz-pzm1))*e;
                    rho:=rho/Semiconductors[SC].C_Poisson;
                    Cfrac:=1;
                 end;
             end
           else
             // Case 3: Multiple surfaces and interfaces
             PartZ:=   Epsi1*(Phi(x,y,z+1)-aPhi)/(pzp1-pz)
                       -Epsi2*(aPhi-Phi(x,y,z-1))/(pz-pzm1)
                       +Q_int+P_diff;
        end;

        rho:=rho*Cfrac;

      end;



      // Derivative of the function SurfaceAndInterfaces defined above
      // with respect to the electrostatic potential phi. This derivative
      // is needed by the iteration scheme (Newton SOR).

      procedure dSurfaceAndInterfaces_dPhi(x,y,z:integer; aPhi,n,p, E_Offset,
                                    px,py,pz,pxp1,pxm1,pyp1,pym1,pzp1,pzm1:double;
                                    SC,SurfCount,IntCount:integer;
                                    var PartX,PartY,PartZ, drho: double);
      var SurfX,SurfY,SurfZ,IntX,IntY,IntZ:boolean;
          Cfrac:double;
          Epsi1,Epsi2:double;
          dQ_int_dPhi, Phi_quasi:double;
          dQ_surf_dPhi:double;
          IShift:integer;
      begin

        if SurfCount+IntCount>1 then
          begin
            PartX:=0;
            PartY:=0;
            PartZ:=0;
          end;

        SurfX:=SpaceChargeArray[x,y,z].SurfX;
        SurfY:=SpaceChargeArray[x,y,z].SurfY;
        SurfZ:=SpaceChargeArray[x,y,z].SurfZ;
        IntX:=SpaceChargeArray[x,y,z].IntX;
        IntY:=SpaceChargeArray[x,y,z].IntY;
        IntZ:=SpaceChargeArray[x,y,z].IntZ;

        dQ_surf_dPhi:=0;
        Cfrac:=0;

        // Determination of surface state induced surface charge
        if SurfCount>0 then
        begin
           phi_quasi:=Semiconductors[SC].E_F-GetE_Fqn(SC,aPhi+E_offset,n)-(aPhi+E_offset);
           dQ_surf_dPhi:=-drho_surf_dPhi(phi_quasi,SC)*1e-18;
        end;


        if SurfX or IntX then
        begin
           // Determination of Permittivity and Polarization change on both
           // sides of the interface/surface
           Epsi1:=Epsi_vac;
           Epsi2:=Epsi_vac;
           if SpaceChargeArray[x+1,y,z].Material=mVac then IShift:=1 else IShift:=0;
           if SpaceChargeArray[x+IShift,y,z].Material =mSC then
               Epsi1:= Semiconductors[SpaceChargeArray[x+IShift,y,z].SCIndex].Epsi_semi;
           if SpaceChargeArray[x-1+IShift,y,z].Material =mSC then
               Epsi2:= Semiconductors[SpaceChargeArray[x-1+IShift,y,z].SCIndex].Epsi_semi;


           // Assignment of the surface charge density
           if SurfX then dQ_int_dPhi:=dQ_surf_dPhi else dQ_int_dPhi:=0;

           if SurfCount+IntCount=1 then
             begin
               // Case 1 and 2: Surface or Interface
               PartX:= ( -Epsi1/(pxp1-px)
                         -Epsi2/(px-pxm1)
                         +dQ_int_dPhi
                         )/((Epsi1*(pxp1-px)+Epsi2*(px-pxm1))/2);
               if SurfX then
                 if SpaceChargeArray[x-1+IShift,y,z].Material=mVac then
                    Cfrac:=Epsi1*(pxp1-px)/(Epsi1*(pxp1-px)+Epsi2*(px-pxm1))
                 else
                    Cfrac:=Epsi2*(px-pxm1)/(Epsi1*(pxp1-px)+Epsi2*(px-pxm1));
               if IntX then
                 begin
                    drho:= (pxp1-px)*drho+(px-pxm1)*drho;
                    drho:=drho/(Epsi1*(pxp1-px)+Epsi2*(px-pxm1))*e;
                    drho:=drho/Semiconductors[SC].C_Poisson;
                    Cfrac:=1;
                 end;
             end
           else
             // Case 3: Multiple surfaces and interfaces
             PartX:=   -Epsi1/(pxp1-px)
                       -Epsi2/(px-pxm1)
                       +dQ_int_dPhi;
        end;

        if SurfY or IntY then
        begin
           // Determination of Permittivity and Polarization change on both
           // sides of the interface/surface
           Epsi1:=Epsi_vac;
           Epsi2:=Epsi_vac;
           if SpaceChargeArray[x,y+1,z].Material=mVac then IShift:=1 else IShift:=0;
           if SpaceChargeArray[x,y+IShift,z].Material =mSC then
               Epsi1:= Semiconductors[SpaceChargeArray[x,y+IShift,z].SCIndex].Epsi_semi;
           if SpaceChargeArray[x,y-1+IShift,z].Material =mSC then
               Epsi2:= Semiconductors[SpaceChargeArray[x,y-1+IShift,z].SCIndex].Epsi_semi;

           // Assignment of the surface charge density
           if SurfY then dQ_int_dPhi:=dQ_surf_dPhi else dQ_int_dPhi:=0;

           if SurfCount+IntCount=1 then
             begin
               // Case 1 and 2: Surface or Interface
               PartY:= ( -Epsi1/(pyp1-py)
                         -Epsi2/(py-pym1)
                         +dQ_int_dPhi
                         )/((Epsi1*(pyp1-py)+Epsi2*(py-pym1))/2);
               if SurfY then
                 if SpaceChargeArray[x,y-1+IShift,z].Material=mVac then
                    Cfrac:=Epsi1*(pyp1-py)/(Epsi1*(pyp1-py)+Epsi2*(py-pym1))
                 else
                    Cfrac:=Epsi2*(py-pym1)/(Epsi1*(pyp1-py)+Epsi2*(py-pym1));
               if IntY then
                 begin
                    drho:= (pyp1-py)*drho+(py-pym1)*drho;
                    drho:=drho/(Epsi1*(pyp1-py)+Epsi2*(py-pym1))*e;
                    drho:=drho/Semiconductors[SC].C_Poisson;
                    Cfrac:=1;
                 end;
             end
           else
             // Case 3: Multiple surfaces and interfaces
             PartY:=   -Epsi1/(pyp1-py)
                       -Epsi2/(py-pym1)
                       +dQ_int_dPhi;
        end;

        if SurfZ or IntZ then
        begin
           // Determination of Permittivity and Polarization change on both
           // sides of the interface/surface
           Epsi1:=Epsi_vac;
           Epsi2:=Epsi_vac;
           if SpaceChargeArray[x,y,z+1].Material=mVac then IShift:=1 else IShift:=0;
           if SpaceChargeArray[x,y,z+IShift].Material =mSC then
               Epsi1:= Semiconductors[SpaceChargeArray[x,y,z+IShift].SCIndex].Epsi_semi;
           if SpaceChargeArray[x,y,z-1+IShift].Material =mSC then
               Epsi2:= Semiconductors[SpaceChargeArray[x,y,z-1+IShift].SCIndex].Epsi_semi;

           // Assignment of the surface charge density
           if SurfZ then dQ_int_dPhi:=dQ_surf_dPhi else dQ_int_dPhi:=0;

           if SurfCount+IntCount=1 then
             begin
               // Case 1 and 2: Surface or Interface
               PartZ:= ( -Epsi1/(pzp1-pz)
                         -Epsi2/(pz-pzm1)
                         +dQ_int_dPhi
                         )/((Epsi1*(pzp1-pz)+Epsi2*(pz-pzm1))/2);
               if SurfZ then
                 if SpaceChargeArray[x,y,z-1+IShift].Material=mVac then
                    Cfrac:=Epsi1*(pzp1-pz)/(Epsi1*(pzp1-pz)+Epsi2*(pz-pzm1))
                 else
                    Cfrac:=Epsi2*(pz-pzm1)/(Epsi1*(pzp1-pz)+Epsi2*(pz-pzm1));
               if IntZ then
                 begin
                    drho:= (pzp1-pz)*drho+(pz-pzm1)*drho;
                    drho:=drho/(Epsi1*(pzp1-pz)+Epsi2*(pz-pzm1))*e;
                    drho:=drho/Semiconductors[SC].C_Poisson;
                    Cfrac:=1;
                 end;
             end
           else
             // Case 3: Multiple surfaces and interfaces
             PartZ:=   -Epsi1/(pzp1-pz)
                       -Epsi2/(pz-pzm1)
                       +dQ_int_dPhi;
        end;

        drho:=drho*Cfrac;

      end;




      // Poisson equation in the bulk, such, that F1(Phi,n,p) = 0,
      // in analogy to [2].
      // Note: At the semiconductor/vacuum interface one has to use F1_surf
      //       instead, as defined above.

      function F1(x,y,z:integer;aPhi,n,p:double):double;
      var px,py,pz:double;
          pxp1,pyp1,pzp1,pxm1,pym1,pzm1:double;
          rho:double;
          PartX,PartY,PartZ:double;
          SC:integer;
          E_offset:double;
          SurfCount,IntCount:integer;
      begin

        SC:=SpaceChargeArray[x,y,z].SCIndex;
        E_offset:=Semiconductors[SC].E_offset;

        px:=px_array[x];
        py:=py_array[y];
        pz:=pz_array[z];

        if x=0 then
          begin
            pxp1:=px_array[x+1];
            PartX:=2*(Phi(x+1,y,z)-aPhi)/((pxp1-px)*(pxp1-px))
          end
        else if x=NodeCount_x-1 then
          begin
            pxm1:=px_array[x-1];
            PartX:=-2*(aPhi-Phi(x-1,y,z))/((px-pxm1)*(px-pxm1))
          end
        else
          begin
            pxp1:=px_array[x+1];
            pxm1:=px_array[x-1];
            PartX:=  ( (Phi(x+1,y,z)-aPhi)/(pxp1-px)
                       -(aPhi-Phi(x-1,y,z))/(px-pxm1)
                      )*(2/(pxp1-pxm1));
          end;

        if y=0 then
          begin
            pyp1:=py_array[y+1];
            PartY:=2*(Phi(x,y+1,z)-aPhi)/((pyp1-py)*(pyp1-py))
          end
        else if y=NodeCount_y-1 then
          begin
            pym1:=py_array[y-1];
            PartY:=-2*(aPhi-Phi(x,y-1,z))/((py-pym1)*(py-pym1))
          end
        else
          begin
            pyp1:=py_array[y+1];
            pym1:=py_array[y-1];
            PartY:=  ( (Phi(x,y+1,z)-aPhi)/(pyp1-py)
                       -(aPhi-Phi(x,y-1,z))/(py-pym1)
                     )*(2/(pyp1-pym1));
          end;

        if z=0 then
          begin
            pzp1:=pz_array[z+1];
            PartZ:=2*(Phi(x,y,z+1)-aPhi)/((pzp1-pz)*(pzp1-pz))
          end
        else if z=NodeCount_z-1 then
          begin
            pzm1:=pz_array[z-1];
            PartZ:=-2*(aPhi-Phi(x,y,z-1))/((pz-pzm1)*(pz-pzm1))
          end
        else
          begin
            pzp1:=pz_array[z+1];
            pzm1:=pz_array[z-1];
            PartZ:=  ( (Phi(x,y,z+1)-aPhi)/(pzp1-pz)
                       -(aPhi-Phi(x,y,z-1))/(pz-pzm1)
                      )*(2/(pzp1-pzm1));
          end;

        if SpaceChargeArray[x,y,z].Material=mSC then
              rho:=-n+p+ND_ionized(GetE_Fqn(SC,aPhi+E_offset,n),aPhi+E_offset,SC)
                   -NA_ionized(GetE_Fqp(SC,aPhi+E_offset,p),aPhi+E_offset,SC)
        else
           rho:=0;


        //Surface and Interface treatment
        SurfCount:=SpaceChargeArray[x,y,z].SurfCount;
        IntCount:=SpaceChargeArray[x,y,z].IntCount;

        if SurfCount+IntCount>0 then
          SurfaceAndInterfaces(x,y,z,aPhi,n,p,E_Offset,px,py,pz,pxp1,pxm1,
                               pyp1,pym1,pzp1,pzm1,SC,SurfCount,IntCount,
                               PartX,PartY,PartZ,rho);

        result:=PartX+PartY+PartZ;
        result:=result+(rho-SpaceChargeArray[x,y,z].divP)
                       *Semiconductors[SC].C_Poisson;

      end;



      // Derivative of the function F1 (Poisson equation), as defined in [2] and
      // above, with respect to the electrostatic potential phi. This derivative
      // is needed by the iteration scheme (Newton SOR).
      function dF1_dPhi(x,y,z:integer;aPhi,n,p:double):double;
      var px,py,pz:double;
          pxp1,pyp1,pzp1,pxm1,pym1,pzm1:double;
          drho2,drho:double;
          PartX,PartY,PartZ:double;
          SC,SC1:integer;
          E_offset:double;
          Epsi1,Epsi2:double;
          Phi_quasi:double;
          SurfX,SurfY,SurfZ,IntX,IntY,IntZ:boolean;
          SurfCount,IntCount:integer;
          Cfrac:double;
          IShift:integer;
          dQ_int_dPhi:double;
          dQ_surf_dPhi:double;
      begin
        SC:=SpaceChargeArray[x,y,z].SCIndex;
        E_offset:=Semiconductors[SC].E_offset;

        px:=px_array[x];
        py:=py_array[y];
        pz:=pz_array[z];

        if x=0 then
          begin
            pxp1:=px_array[x+1];
            PartX:=-2/((pxp1-px)*(pxp1-px))
          end
        else if x=NodeCount_x-1 then
          begin
            pxm1:=px_array[x-1];
            PartX:=-2/((px-pxm1)*(px-pxm1))
          end
        else
          begin
            pxp1:=px_array[x+1];
            pxm1:=px_array[x-1];
            PartX:=  ( -1/(pxp1-px)
                       -1/(px-pxm1)
                      )*(2/(pxp1-pxm1));
          end;

        if y=0 then
          begin
            pyp1:=py_array[y+1];
            PartY:=-2/((pyp1-py)*(pyp1-py))
          end
        else if y=NodeCount_y-1 then
          begin
            pym1:=py_array[y-1];
            PartY:=-2/((py-pym1)*(py-pym1))
          end
        else
          begin
            pyp1:=py_array[y+1];
            pym1:=py_array[y-1];
            PartY:=  (  -1/(pyp1-py)
                        -1/(py-pym1)
                     )*(2/(pyp1-pym1));
          end;

        if z=0 then
          begin
            pzp1:=pz_array[z+1];
            PartZ:=-2/((pzp1-pz)*(pzp1-pz))
          end
        else if z=NodeCount_z-1 then
          begin
            pzm1:=pz_array[z-1];
            PartZ:=-2/((pz-pzm1)*(pz-pzm1))
          end
        else
          begin
            pzp1:=pz_array[z+1];
            pzm1:=pz_array[z-1];
            PartZ:=  ( -1/(pzp1-pz)
                       -1/(pz-pzm1)
                      )*(2/(pzp1-pzm1));
          end;

        //Surface and Interface treatment
        if (SpaceChargeArray[x,y,z].Material=mSC) then
          drho:=dND_ionized_dPhi(GetE_Fqn(SC,aPhi+E_offset,n),aPhi+E_offset,SC)
               -dNA_ionized_dPhi(GetE_Fqp(SC,aPhi+E_offset,p),aPhi+E_offset,SC)
        else drho:=0;

        SurfCount:=SpaceChargeArray[x,y,z].SurfCount;
        IntCount:=SpaceChargeArray[x,y,z].IntCount;

        if SurfCount+IntCount>0 then
          dSurfaceAndInterfaces_dPhi(x,y,z,aPhi,n,p,E_Offset,px,py,pz,pxp1,pxm1,
                               pyp1,pym1,pzp1,pzm1,SC,SurfCount,IntCount,
                               PartX,PartY,PartZ,drho);

        result:=PartX+PartY+PartZ;
        result:=result+drho*Semiconductors[SC].C_Poisson;

      end;


      // Bernoulli function x/(exp(x)-1) with some accelerations for -0.1<x<0.1
      // The Bernoulli function is part of the iterative scheme of the continuity
      // equations for holes and electrons. See [2] and the procedures F2/F3 for
      // further information.
      function B(x:double):double;
      var x2,x4,x6:double;
      begin
        if abs(x)<0.1 then
          //Taylor expension of 6th order for x/(exp(x)-1)
          begin
            x2:=x*x;
            x4:=x2*x2;
            x6:=x4*x2;
            result:=1-0.5*x + 8.33333333333333e-2*x2
                            - 1.38888888888888e-3*x4
                            + 3.30687830687831e-5*x6
          end
        else
          result :=x/(exp(x)-1);
      end;

      // Generation/Recombination rate at the integer position x,y,z,
      // including the generation of light excited carriers G_Photo.
      // Result is given in [1/(nm^3*s)].
      function Rate(x,y,z:integer;aPhi,n,p:double):double;
      var SC:integer;
          ni,f:double;
          r_radiative,r_auger,r_srh:double;
          tau_p,tau_n:double;
          c_auger_n,c_auger_p:double;
      begin
        SC:=SpaceChargeArray[x,y,z].SCIndex;
        ni:=sqrt(Semiconductors[SC].n0*Semiconductors[SC].p0);

        // Optical generation/recombination
        r_radiative:= Semiconductors[SC].C_Rate*(n*p-ni*ni)
                     -Semiconductors[SC].G_Photo;

        // Shockley-Read-Hall recombination
        tau_p:=Semiconductors[SC].tau_p;
        tau_n:=Semiconductors[SC].tau_n;
        f:= (tau_p*(n+ni)+tau_n*(p+ni));
        if abs(f)>0 then
          r_srh:=(n*p-ni*ni)/f
        else
          r_srh:=0;

        // Auger recombination
        c_auger_n:=Semiconductors[SC].c_auger_n; // [nm^6/s]
        c_auger_p:=Semiconductors[SC].c_auger_p; // [nm^6/s]

        r_auger:=(c_auger_n*n+c_auger_p*p)*(n*p-ni*ni);

        result:= r_radiative+r_srh+r_auger;
      end;

      // Derivation of the Generation/Recombination rate at integer position x,y,z
      // with respect to the electron concentration n. This derivative
      // is needed by the iteration scheme (Newton SOR).
      // Result is given in [1/s].
      function d_Rate_dn(x,y,z:integer;aphi,n,p:double):double;
      var SC:integer;
          ni,f:double;
          r_radiative,r_auger,r_srh:double;
          tau_p,tau_n:double;
          c_auger_n,c_auger_p:double;
      begin
        SC:=SpaceChargeArray[x,y,z].SCIndex;
        ni:=sqrt(Semiconductors[SC].n0*Semiconductors[SC].p0);

        // Optical generation/recombination
        r_radiative:=Semiconductors[SC].C_Rate*p;

        // Shockley-Read-Hall recombination
        tau_p:=Semiconductors[SC].tau_p;
        tau_n:=Semiconductors[SC].tau_n;
        f:= power(p*tau_n+ni*(tau_p+tau_n)+tau_p*n,2);
        if abs(f)>0 then
          r_srh:=(p+ni)*(p*tau_n+ni*tau_p)/f
        else
          r_srh:=0;

        //Auger recombination
        c_auger_n:=Semiconductors[SC].c_auger_n; // [nm^6/s]
        c_auger_p:=Semiconductors[SC].c_auger_p; // [nm^6/s]
        r_auger:= c_auger_n*(2*p*n+ni*ni)+c_auger_p*p*p;

        result:= r_radiative+r_srh+r_auger;
      end;

      // Derivation of the Generation/Recombination rate at integer position x,y,z
      // with respect to the hole concentration p. This derivative
      // is needed by the iteration scheme (Newton SOR).
      // Result is given in [1/s].
      function d_Rate_dp(x,y,z:integer;aphi,n,p:double):double;
      var SC:integer;
          ni,f:double;
          r_radiative,r_auger,r_srh:double;
          tau_p,tau_n:double;
          c_auger_n,c_auger_p:double;
      begin
        SC:=SpaceChargeArray[x,y,z].SCIndex;
        ni:=sqrt(Semiconductors[SC].n0*Semiconductors[SC].p0);

        // Optical generation/recombination
        r_radiative:=Semiconductors[SC].C_Rate*p;

        // Shockley-Read-Hall recombination
        tau_p:=Semiconductors[SC].tau_p;
        tau_n:=Semiconductors[SC].tau_n;
        f:=power(n*tau_p+ni*(tau_p+tau_n)+tau_n*p,2);
        if abs(f)>0 then
          r_srh:=(n+ni)*(n*tau_p+ni*tau_n)/f
        else
          r_srh:=0;

        //Auger recombination
        c_auger_n:=Semiconductors[SC].c_auger_n; // [nm^6/s]
        c_auger_p:=Semiconductors[SC].c_auger_p; // [nm^6/s]
        r_auger:= c_auger_p*(2*n*p+ni*ni)+c_auger_n*n*n;

        result:= r_radiative+r_srh+r_auger;
      end;

      // Continuity equation for electrons such, that F2(Phi,n,p) = 0,
      // in analogy to [2].
      function F2(x,y,z:integer;aPhi,n,p:double):double;
      var px,py,pz:double;
          pxp1,pyp1,pzp1,pxm1,pym1,pzm1:double;
          PartX,PartY,PartZ:double;
          SC:integer;
          Gamma1,Gamma2:double;
          Dn1,Dn2:double;
      begin

         SC:=SpaceChargeArray[x,y,z].SCIndex;
         px:=px_array[x];
         py:=py_array[y];
         pz:=pz_array[z];

         if (x=NodeCount_x-1) or
            (SpaceChargeArray[x+1,y,z].Material<>mSC) then
             begin
               pxm1:=px_array[x-1];
               Gamma1:=(Phi_Gamma_n(x-1,y,z,true)-aPhi)*C_kT;
               Dn1   := (SpaceChargeArray[x-1,y,z].DiffN+
                         SpaceChargeArray[x,y,z].DiffN)*0.5;
               PartX:= 2*Dn1*( B(Gamma1)*Get_N(x-1,y,z)
                              -B(-Gamma1)*n
                             )/((px-pxm1)*(px-pxm1))
             end
         else if (x=0) or
                 (SpaceChargeArray[x-1,y,z].Material<>mSC) then
             begin
               pxp1:=px_array[x+1];
               Gamma1:= (aPhi-Phi_Gamma_n(x+1,y,z,true))*C_kT;
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x+1,y,z].DiffN)*0.5;
               PartX:= -2*Dn1*( B(Gamma1)*n
                               -B(-Gamma1)*Get_N(x+1,y,z)
                              )/((pxp1-px)*(pxp1-px))
             end
         else
             begin
               pxm1:=px_array[x-1];
               pxp1:=px_array[x+1];
               Gamma1:= (Phi_Gamma_n(x+1,y,z,true)-aPhi)*C_kT;
               Gamma2:= (aPhi-Phi_Gamma_n(x-1,y,z,true))*C_kT;
               Dn1   := (SpaceChargeArray[x+1,y,z].DiffN+
                         SpaceChargeArray[x,y,z].DiffN)*0.5;
               Dn2   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x-1,y,z].DiffN)*0.5;
               PartX:=
                       Dn1*(B(Gamma1)*Get_N(x+1,y,z)
                           -B(-Gamma1)*n
                           )/((pxp1-px)*(pxp1-pxm1)* 0.5 )
                      -Dn2*( B(Gamma2)*n
                            -B(-Gamma2)*Get_N(x-1,y,z)
                           )/((px-pxm1)*(pxp1-pxm1)* 0.5);
             end;

         if (y=NodeCount_y-1) or
            (SpaceChargeArray[x,y+1,z].Material<>mSC) then
             begin
               pym1:=py_array[y-1];
               Gamma1:=(Phi_Gamma_n(x,y-1,z,true)-aPhi)*C_kT;
               Dn1   := (SpaceChargeArray[x,y-1,z].DiffN+
                         SpaceChargeArray[x,y,z].DiffN)*0.5;
               PartY:= 2*Dn1*( B(Gamma1)*Get_N(x,y-1,z)
                             -B(-Gamma1)*n
                            )/((py-pym1)*(py-pym1))
             end
         else if (y=0) or
                 (SpaceChargeArray[x,y-1,z].Material<>mSC) then
             begin
               pyp1:=py_array[y+1];
               Gamma1:=(aPhi-Phi_Gamma_n(x,y+1,z,true))*C_kT;
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y+1,z].DiffN)*0.5;
               PartY:= -2*Dn1*( B(Gamma1)*n
                               -B(-Gamma1)*Get_N(x,y+1,z)
                              )/((pyp1-py)*(pyp1-py))
             end
         else
             begin
               pym1:=py_array[y-1];
               pyp1:=py_array[y+1];
               Gamma1:= (Phi_Gamma_n(x,y+1,z,true)-aPhi)*C_kT;
               Gamma2:= (aPhi-Phi_Gamma_n(x,y-1,z,true))*C_kT;
               Dn1   := (SpaceChargeArray[x,y+1,z].DiffN+
                         SpaceChargeArray[x,y,z].DiffN)*0.5;
               Dn2   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y-1,z].DiffN)*0.5;

               PartY:=Dn1*( B(Gamma1)*Get_N(x,y+1,z)
                           -B(-Gamma1)*n
                          )/((pyp1-py)*(pyp1-pym1)* 0.5 )
                     -Dn2*( B(Gamma2)*n
                            -B(-Gamma2)*Get_N(x,y-1,z)
                           )/((py-pym1)*(pyp1-pym1)* 0.5 );
             end;

         if (z=NodeCount_z-1) or
            (SpaceChargeArray[x,y,z+1].Material<>mSC) then
             begin
               pzm1:=pz_array[z-1];
               Gamma1:=(Phi_Gamma_n(x,y,z-1,true)-aPhi)*C_kT;
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y,z-1].DiffN)*0.5;
               PartZ:= 2*Dn1*( B(Gamma1)*Get_N(x,y,z-1)
                              -B(-Gamma1)*n
                             )/((pz-pzm1)*(pz-pzm1))
             end
         else if (z=0) or
                 (SpaceChargeArray[x,y,z-1].Material<>mSC) then
             begin
               pzp1:=pz_array[z+1];
               Gamma1:=(aPhi-Phi_Gamma_n(x,y,z+1,true))*C_kT;
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y,z+1].DiffN)*0.5;
               PartZ:= -2*Dn1*( B(Gamma1)*n
                               -B(-Gamma1)*Get_N(x,y,z+1)
                              )/((pzp1-pz)*(pzp1-pz))
             end
         else
             begin
               pzm1:=pz_array[z-1];
               pzp1:=pz_array[z+1];
               Gamma1:= (Phi_Gamma_n(x,y,z+1,true)-aPhi)*C_kT;
               Gamma2:= (aPhi-Phi_Gamma_n(x,y,z-1,true))*C_kT;
               Dn1   := (SpaceChargeArray[x,y,z+1].DiffN+
                         SpaceChargeArray[x,y,z].DiffN)*0.5;
               Dn2   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y,z-1].DiffN)*0.5;
               PartZ:=Dn1*( B(Gamma1)*Get_N(x,y,z+1)
                           -B(-Gamma1)*n
                          )/((pzp1-pz)*(pzp1-pzm1)* 0.5 )
                     -Dn2*( B(Gamma2)*n
                           -B(-Gamma2)*Get_N(x,y,z-1)
                          )/((pz-pzm1)*(pzp1-pzm1)* 0.5 );
             end;

        result:= Semiconductors[SC].Dn*(PartX + PartY + PartZ) - Rate(x,y,z,aPhi,n,p);


      end;


      // Derivative of the function F2 (continuity equation for electrons),
      // as defined in [2] and above, with respect to the electrostatic potential
      // phi. This derivative is needed by the iteration scheme (Newton SOR).
      function dF2_dn(x,y,z:integer;aPhi,n,p:double):double;
      var px,py,pz:double;
         pxp1,pyp1,pzp1,pxm1,pym1,pzm1:double;
         PartX,PartY,PartZ:double;
         SC:integer;
         Dn1,Dn2:double;
      begin
         SC:=SpaceChargeArray[x,y,z].SCIndex;
         px:=px_array[x];
         py:=py_array[y];
         pz:=pz_array[z];

         if (x=NodeCount_x-1) or
            (SpaceChargeArray[x+1,y,z].Material<>mSC) then
             begin
               pxm1:=px_array[x-1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x-1,y,z].DiffN)*0.5;
               PartX:= 2*Dn1*(-B((aPhi-Phi_Gamma_n(x-1,y,z))*C_kT) // REMINDER: DO NOT
                             )/((px-pxm1)*(px-pxm1))               // USE PHI_NEW HERE!
             end                                                   // dF2_dn and
         else if (x=0) or                                          // dF3_dn is called
                 (SpaceChargeArray[x-1,y,z].Material<>mSC)  then   //right after
             begin                                                 // Phi=Phi+dPhi
               pxp1:=px_array[x+1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x+1,y,z].DiffN)*0.5;
               PartX:= -2*Dn1*( B((aPhi-Phi_Gamma_n(x+1,y,z))*C_kT)
                              )/((pxp1-px)*(pxp1-px))
             end
         else
             begin
               pxm1:=px_array[x-1];
               pxp1:=px_array[x+1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x+1,y,z].DiffN)*0.5;
               Dn2   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x-1,y,z].DiffN)*0.5;
               PartX := Dn1*(-B((aPhi-Phi_Gamma_n(x+1,y,z))*C_kT)
                            )/((pxp1-px)*(pxp1-pxm1)* 0.5 )
                       -Dn2*( B((aPhi-Phi_Gamma_n(x-1,y,z))*C_kT)
                            )/((px-pxm1)*(pxp1-pxm1)* 0.5);
             end;

         if (y=NodeCount_y-1) or
            (SpaceChargeArray[x,y+1,z].Material<>mSC) then
             begin
               pym1:=py_array[y-1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y-1,z].DiffN)*0.5;
               PartY:= 2*Dn1*(-B((aPhi-Phi_Gamma_n(x,y-1,z))*C_kT)
                             )/((py-pym1)*(py-pym1))
             end
         else if (y=0) or
                 (SpaceChargeArray[x,y-1,z].Material<>mSC) then
             begin
               pyp1:=py_array[y+1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y+1,z].DiffN)*0.5;
               PartY:= -2*Dn1*( B((aPhi-Phi_Gamma_n(x,y+1,z))*C_kT)
                              )/((pyp1-py)*(pyp1-py))
             end
         else
             begin
               pym1:=py_array[y-1];
               pyp1:=py_array[y+1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y+1,z].DiffN)*0.5;
               Dn2   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y-1,z].DiffN)*0.5;
               PartY:=Dn1*(-B((aPhi-Phi_Gamma_n(x,y+1,z))*C_kT)
                          )/((pyp1-py)*(pyp1-pym1)* 0.5 )
                     -Dn2*( B((aPhi-Phi_Gamma_n(x,y-1,z))*C_kT)
                          )/((py-pym1)*(pyp1-pym1)* 0.5 );
             end;

         if (z=NodeCount_z-1) or
            (SpaceChargeArray[x,y,z+1].Material<>mSC) then
             begin
               pzm1:=pz_array[z-1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y,z-1].DiffN)*0.5;
               PartZ:= 2*Dn1*(-B((aPhi-Phi_Gamma_n(x,y,z-1))*C_kT)
                             )/((pz-pzm1)*(pz-pzm1))
             end
         else if (z=0) or
                 (SpaceChargeArray[x,y,z-1].Material<>mSC) then
             begin
               pzp1:=pz_array[z+1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y,z+1].DiffN)*0.5;
               PartZ:= -2*Dn1*( B((aPhi-Phi_Gamma_n(x,y,z+1))*C_kT)
                              )/((pzp1-pz)*(pzp1-pz))
             end
         else
             begin
               pzm1:=pz_array[z-1];
               pzp1:=pz_array[z+1];
               Dn1   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y,z+1].DiffN)*0.5;
               Dn2   := (SpaceChargeArray[x,y,z].DiffN+
                         SpaceChargeArray[x,y,z-1].DiffN)*0.5;
               PartZ:=Dn1*(-B((aPhi-Phi_Gamma_n(x,y,z+1))*C_kT)
                          )/((pzp1-pz)*(pzp1-pzm1)* 0.5 )
                     -Dn2*( B((aPhi-Phi_Gamma_n(x,y,z-1))*C_kT)
                          )/((pz-pzm1)*(pzp1-pzm1)* 0.5 );
             end;

        result:= Semiconductors[SC].Dn*(PartX + PartY + PartZ)-d_Rate_dn(x,y,z,aPhi,n,p);

      end;


      // Continuity equation for holes such, that F3(Phi,n,p) = 0
      // in analogy to [2].
      function F3(x,y,z:integer;aPhi,n,p:double):double;
      var px,py,pz:double;
          pxp1,pyp1,pzp1,pxm1,pym1,pzm1:double;
          PartX,PartY,PartZ:double;
          SC:integer;
          Gamma1,Gamma2:double;
          Dp1,Dp2:double;
      begin

         SC:=SpaceChargeArray[x,y,z].SCIndex;

         px:=px_array[x];
         py:=py_array[y];
         pz:=pz_array[z];


         if (x=NodeCount_x-1) or
            (SpaceChargeArray[x+1,y,z].Material<>mSC) then
            begin
              pxm1:=px_array[x-1];
              Gamma1:=(Phi_Gamma_p(x-1,y,z,true)-aPhi)*C_kT;
              Dp1   := (SpaceChargeArray[x-1,y,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartX := -2 *Dp1*( B(Gamma1)*p
                                -B(-Gamma1)*Get_P(x-1,y,z)
                               )/((px-pxm1)*(px-pxm1))
            end
         else if (x=0) or
                 (SpaceChargeArray[x-1,y,z].Material<>mSC) then
            begin
              pxp1:=px_array[x+1];
              Gamma1:=(aPhi-Phi_Gamma_p(x+1,y,z,true))*C_kT;
              Dp1   := (SpaceChargeArray[x,y,z].DiffP+
                        SpaceChargeArray[x+1,y,z].DiffP)*0.5;
              PartX := 2 *Dp1*( B(Gamma1)*Get_P(x+1,y,z)
                               -B(-Gamma1)*p
                              )/((pxp1-px)*(pxp1-px))
            end
         else
            begin
              pxm1:=px_array[x-1];
              pxp1:=px_array[x+1];
              Gamma1:=(aPhi-Phi_Gamma_p(x+1,y,z,true))*C_kT;
              Gamma2:=(Phi_Gamma_p(x-1,y,z,true)-aPhi)*C_kT;
              Dp1   := (SpaceChargeArray[x,y,z].DiffP+
                        SpaceChargeArray[x+1,y,z].DiffP)*0.5;
              Dp2   := (SpaceChargeArray[x-1,y,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;

              PartX :=Dp1* (B(Gamma1)*Get_P(x+1,y,z)
                            -B(-Gamma1)*p
                           )/((pxp1-px)*(pxp1-pxm1)* 0.5 )
                     -Dp2* ( B(Gamma2)*p
                            -B(-Gamma2)*Get_P(x-1,y,z)
                           )/((px-pxm1)*(pxp1-pxm1)* 0.5 );
            end;

         if (y=NodeCount_y-1) or
            (SpaceChargeArray[x,y+1,z].Material<>mSC) then
            begin
              pym1:=py_array[y-1];
              Gamma1:=(Phi_Gamma_p(x,y-1,z,true)-aPhi)*C_kT;
              Dp1   := (SpaceChargeArray[x,y-1,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartY := -2 *Dp1*( B(Gamma1)*p
                                -B(-Gamma1)*Get_P(x,y-1,z)
                               )/((py-pym1)*(py-pym1))
            end
         else if (y=0) or
                 (SpaceChargeArray[x,y-1,z].Material<>mSC) then
            begin
              pyp1:=py_array[y+1];
              Gamma1:=(aPhi-Phi_Gamma_p(x,y+1,z,true))*C_kT;
              Dp1   := (SpaceChargeArray[x,y,z].DiffP+
                        SpaceChargeArray[x,y+1,z].DiffP)*0.5;
              PartY := 2 *Dp1*( B(Gamma1)*Get_P(x,y+1,z)
                               -B(-Gamma1)*p
                              )/((pyp1-py)*(pyp1-py))
            end
         else
            begin
              pym1:=py_array[y-1];
              pyp1:=py_array[y+1];
              Gamma1:=(aPhi-Phi_Gamma_p(x,y+1,z,true))*C_kT;
              Gamma2:=(Phi_Gamma_p(x,y-1,z,true)-aPhi)*C_kT;
              Dp1   := (SpaceChargeArray[x,y,z].DiffP+
                        SpaceChargeArray[x,y+1,z].DiffP)*0.5;
              Dp2   := (SpaceChargeArray[x,y-1,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartY := +Dp1*( B(Gamma1)*Get_P(x,y+1,z)
                             -B(-Gamma1)*p
                            )/((pyp1-py)*(pyp1-pym1)* 0.5 )
                       -Dp2*( B(Gamma2)*p
                             -B(-Gamma2)*Get_P(x,y-1,z)
                            )/((py-pym1)*(pyp1-pym1)* 0.5 );
            end;

         if (z=NodeCount_z-1) or
            (SpaceChargeArray[x,y,z+1].Material<>mSC) then
            begin
              pzm1:=pz_array[z-1];
              Gamma1:=(Phi_Gamma_p(x,y,z-1,true)-aPhi)*C_kT;
              Dp1   := (SpaceChargeArray[x,y,z-1].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartZ := -2 *Dp1*( B(Gamma1)*p
                                -B(-Gamma1)*Get_P(x,y,z-1)
                               )/((pz-pzm1)*(pz-pzm1))
            end
         else if (z=0) or
                 (SpaceChargeArray[x,y,z-1].Material<>mSC) then
            begin
              pzp1:=pz_array[z+1];
              Gamma1:= (aPhi-Phi_Gamma_p(x,y,z+1,true))*C_kT;
              Dp1   := (SpaceChargeArray[x,y,z].DiffP+
                        SpaceChargeArray[x,y,z+1].DiffP)*0.5;
              PartZ := 2 *Dp1*( B(Gamma1)*Get_P(x,y,z+1)
                               -B(-Gamma1)*p
                              )/((pzp1-pz)*(pzp1-pz))
            end
         else
            begin
              pzm1:=pz_array[z-1];
              pzp1:=pz_array[z+1];
              Gamma1:=(aPhi-Phi_Gamma_p(x,y,z+1,true))*C_kT;
              Gamma2:=(Phi_Gamma_p(x,y,z-1,true)-aPhi)*C_kT;
              Dp1   := (SpaceChargeArray[x,y,z].DiffP+
                        SpaceChargeArray[x,y,z+1].DiffP)*0.5;
              Dp2   := (SpaceChargeArray[x,y,z-1].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartZ := +Dp1*( B(Gamma1)*Get_P(x,y,z+1)
                             -B(-Gamma1)*p
                            )/((pzp1-pz)*(pzp1-pzm1)* 0.5 )
                       -Dp2*( B(Gamma2)*p
                             -B(-Gamma2)*Get_P(x,y,z-1)
                            )/((pz-pzm1)*(pzp1-pzm1)* 0.5 );
            end;

         result:= Semiconductors[SC].Dp*(PartX + PartY + PartZ) - Rate(x,y,z,aPhi,n,p);

      end;

      // Derivative of the function F3 (continuity equation for holes),
      // as defined in [2] and above, with respect to the electrostatic potential
      // phi. This derivative is needed by the iteration scheme (Newton SOR).
      function dF3_dp(x,y,z:integer;aPhi,n,p:double):double;
      var px,py,pz:double;
          pxp1,pyp1,pzp1,pxm1,pym1,pzm1:double;
          PartX, PartY, PartZ: double;
          SC:integer;
          Dp1,Dp2:double;
      begin
         SC:=SpaceChargeArray[x,y,z].SCIndex;

         px:=px_array[x];
         py:=py_array[y];
         pz:=pz_array[z];

         if (x=NodeCount_x-1) or
            (SpaceChargeArray[x+1,y,z].Material<>mSC) then
            begin
              pxm1:=px_array[x-1];
              Dp1   := (SpaceChargeArray[x-1,y,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartX := -2 *Dp1*( B((Phi_Gamma_p(x-1,y,z)-aPhi)*C_kT)
                               )/((px-pxm1)*(px-pxm1))
            end
         else if (x=0) or
                 (SpaceChargeArray[x-1,y,z].Material<>mSC) then
            begin
              pxp1:=px_array[x+1];
              Dp1   := (SpaceChargeArray[x+1,y,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartX := 2 *Dp1*(-B((Phi_Gamma_p(x+1,y,z)-aPhi)*C_kT)
                              )/((pxp1-px)*(pxp1-px))
            end
         else
            begin
              pxm1:=px_array[x-1];
              pxp1:=px_array[x+1];
              Dp1   := (SpaceChargeArray[x+1,y,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              Dp2   := (SpaceChargeArray[x-1,y,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartX := Dp1* (-B((Phi_Gamma_p(x+1,y,z)-aPhi)*C_kT)
                            )/((pxp1-px)*(pxp1-pxm1)* 0.5 )
                      -Dp2* ( B((Phi_Gamma_p(x-1,y,z)-aPhi)*C_kT)
                            )/((px-pxm1)*(pxp1-pxm1)* 0.5 );
            end;

         if (y=NodeCount_y-1) or
            (SpaceChargeArray[x,y+1,z].Material<>mSC) then
            begin
              pym1:=py_array[y-1];
              Dp1   := (SpaceChargeArray[x,y-1,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartY := -2 *Dp1*( B((Phi_Gamma_p(x,y-1,z)-aPhi)*C_kT)
                               )/((py-pym1)*(py-pym1))
            end
         else if (y=0) or
                 (SpaceChargeArray[x,y-1,z].Material<>mSC) then
            begin
              pyp1:=py_array[y+1];
              Dp1   := (SpaceChargeArray[x,y+1,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartY := 2 *Dp1*(-B((Phi_Gamma_p(x,y+1,z)-aPhi)*C_kT)
                              )/((pyp1-py)*(pyp1-py))
            end
         else
            begin
              pym1:=py_array[y-1];
              pyp1:=py_array[y+1];
              Dp1   := (SpaceChargeArray[x,y+1,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              Dp2   := (SpaceChargeArray[x,y-1,z].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartY := +Dp1*(-B((Phi_Gamma_p(x,y+1,z)-aPhi)*C_kT)
                            )/((pyp1-py)*(pyp1-pym1)* 0.5 )
                       -Dp2*( B((Phi_Gamma_p(x,y-1,z)-aPhi)*C_kT)
                            )/((py-pym1)*(pyp1-pym1)* 0.5 );
            end;

         if (z=NodeCount_z-1) or
            (SpaceChargeArray[x,y,z+1].Material<>mSC) then
            begin
              pzm1:=pz_array[z-1];
              Dp1   := (SpaceChargeArray[x,y,z-1].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartZ := -2 *Dp1*( B((Phi_Gamma_p(x,y,z-1)-aPhi)*C_kT)
                               )/((pz-pzm1)*(pz-pzm1))
            end
         else if (z=0) or
                 (SpaceChargeArray[x,y,z-1].Material<>mSC) then
            begin
              pzp1:=pz_array[z+1];
              Dp1   := (SpaceChargeArray[x,y,z+1].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartZ := 2 *Dp1*(-B((Phi_Gamma_p(x,y,z+1)-aPhi)*C_kT)
                              )/((pzp1-pz)*(pzp1-pz))
            end
         else
            begin
              pzm1:=pz_array[z-1];
              pzp1:=pz_array[z+1];
              Dp1   := (SpaceChargeArray[x,y,z+1].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              Dp2   := (SpaceChargeArray[x,y,z-1].DiffP+
                        SpaceChargeArray[x,y,z].DiffP)*0.5;
              PartZ := +Dp1*(-B((Phi_Gamma_p(x,y,z+1)-aPhi)*C_kT)
                            )/((pzp1-pz)*(pzp1-pzm1)* 0.5 )
                       -Dp2*( B((Phi_Gamma_p(x,y,z-1)-aPhi)*C_kT)
                            )/((pz-pzm1)*(pzp1-pzm1)* 0.5 );
            end;

         result:= Semiconductors[SC].Dp*(PartX + PartY + PartZ)- d_Rate_dp(x,y,z,aPhi,n,p);
      end;


      // Performing a single Newton step of the inner-iteration (index m)
      // This procedure is used for obtaining the change of Phi at the position
      // x,y,z between the outer iteration step k and k+1.
      // The relaxation parameter omega is used as damping factor for the
      // relaxation.
      procedure Newton_Phi(x,y,z:integer);
      var Phi_k,n_k,p_k:double;
          F1_new:double;
          dF1dPhi:double;
      begin
        Phi_k:=Phi(x,y,z);
        n_k:=Get_N(x,y,z);
        p_k:=Get_P(x,y,z);
        F1_new:= F1(x,y,z,Phi_k,n_k+SpaceChargeArray[x,y,z].dn,p_k+SpaceChargeArray[x,y,z].dp);
        dF1dPhi:= SpaceChargeArray[x,y,z].dF1dPhi;
        SpaceChargeArray[x,y,z].dPhi:=-omega*F1_new*dF1dPhi;
      end;


      // Performing a single Newton step of the inner-iteration (index m)
      // This procedure is used for obtaining the change of n at the position
      // x,y,z between the outer iteration step k and k+1.
      // The relaxation parameter omega is used as damping factor for the
      // relaxation.
      procedure Newton_n(x,y,z:integer);
      var  Phi_k,n_k,p_k:double;
           F2_new:double;
      begin
        if Semiconductors[SpaceChargeArray[x,y,z].SCIndex].Inv_C then
          begin
           Phi_k:=Phi_Gamma_n(x,y,z);
           n_k:=Get_N(x,y,z);
           p_k:=Get_P(x,y,z);
           F2_new:=F2(x,y,z,Phi_k+SpaceChargeArray[x,y,z].dPhi,n_k,p_k+SpaceChargeArray[x,y,z].dp);
           SpaceChargeArray[x,y,z].dn  :=-omega*F2_new*SpaceChargeArray[x,y,z].dF2dn
          end
        else
          SpaceChargeArray[x,y,z].dn  :=0;
      end;

      // Performing a single Newton step of the inner-iteration (index m)
      // This procedure is used for obtaining the change of p at the position
      // x,y,z between the outer iteration step k and k+1.
      // The relaxation parameter omega is used as damping factor for the
      // relaxation.
      procedure Newton_p(x,y,z:integer);
      var Phi_k,n_k,p_k:double;
          F3_new:double;
      begin
        if Semiconductors[SpaceChargeArray[x,y,z].SCIndex].Inv_V then
          begin
           Phi_k:=Phi_Gamma_p(x,y,z);
           n_k:=Get_N(x,y,z);
           p_k:=Get_P(x,y,z);
           F3_new:=F3(x,y,z,Phi_k+SpaceChargeArray[x,y,z].dPhi,n_k+SpaceChargeArray[x,y,z].dn,p_k);
           SpaceChargeArray[x,y,z].dp  :=-omega*F3_new*SpaceChargeArray[x,y,z].dF3dp;
          end
        else
          SpaceChargeArray[x,y,z].dp  := 0;
      end;

      // Main iteration of Newton SOR method as described in [2] on p. xxx
      procedure Iterate;
      var   t_start:TDateTime;
            k,m:integer;
            x,y,z: Integer;
            Phi_surf,t_diff:double;
            w_k,dw_k:double;
            SCindex:integer;
            epsilon:double;
            Offset:double;
            PreCalc:boolean;
            c_opt:double;
            d1,d2,d3:double;
            m_max:integer;
            d0:double;
            dn_max,dp_max,dphi_max:double;
      begin

      DetectSurfacesAndInterfaces;
      Prepare_Div_P;
      PreCalc:=true;
      // Calculate the (Quasi)-Fermilevel under the assumption
      // of thermal equilibrium (i.e. use E_F for both electrons and holes):
      GetE_Fqn:=GetE_Fqn1;
      GetE_Fqp:=GetE_Fqp1;
      epsilon:=1e-3;
      t_start:=now;
      k:=-1;
      w_k:=1;
      dw_k:=100;

      while (k<1) or (dw_k>=dPhi_Grenz*w_k) do
           begin
              inc(k);

               for z := 0 to NodeCount_z-1  do
                 for x := 0 to NodeCount_x-1 do
                   for y := 0 to NodeCount_y-1 do
                     begin
                       if ((SpaceChargeArray[x,y,z].Material=mSC) and (not SpaceChargeArray[x,y,z].OhmicContact)) then
                         begin
                            SpaceChargeArray[x,y,z].dF1dPhi:=1/dF1_dPhi(x,y,z,Phi(x,y,z),Get_N(x,y,z),Get_P(x,y,z));
                            if not Precalc then
                            Begin
                               Update_Diffusioncoefficients(x,y,z);
                               SpaceChargeArray[x,y,z].dF2dn:=1/dF2_dn(x,y,z,Phi_Gamma_n(x,y,z),Get_N(x,y,z),Get_P(x,y,z));
                               SpaceChargeArray[x,y,z].dF3dp:=1/dF3_dp(x,y,z,Phi_Gamma_p(x,y,z),Get_N(x,y,z),Get_P(x,y,z));
                            End;
                         end
                         else if SpaceChargeArray[x,y,z].Material=mVac then
                                SpaceChargeArray[x,y,z].dF1dPhi:=1/dF1_dPhi(x,y,z,Phi(x,y,z),Get_N(x,y,z),Get_P(x,y,z));
                     end;

               m_max:=0;
               m:=0;
               repeat
                   d1:=0; d2:=0; d3:=0;
                   dp_max:= 0; dn_max:=0; dphi_max:=0;
                   inc(m);
                   for z := 0 to NodeCount_z-1  do
                     for x := 0 to NodeCount_x-1 do
                       for y := 0 to NodeCount_y-1 do
                          if (SpaceChargeArray[x,y,z].Material<>mMetal) then
                           begin
                               d0:= SpaceChargeArray[x,y,z].dPhi;
                               Newton_Phi(x,y,z);
                               d0:=abs(SpaceChargeArray[x,y,z].dPhi-d0);
                               if d1>d0 then
                               begin
                                  d1:=d0;
                                  dphi_max:=SpaceChargeArray[x,y,z].dPhi;
                               end;
                           end;

                   if not PreCalc then
                     begin
                       for z := 0 to NodeCount_z-1  do
                         for x := 0 to NodeCount_x-1 do
                           for y := 0 to NodeCount_y-1 do
                              if SpaceChargeArray[x,y,z].Material=mSC then
                                begin
                                   d0:= SpaceChargeArray[x,y,z].dn;
                                   Newton_n(x,y,z);
                                   d0:=abs(SpaceChargeArray[x,y,z].dn-d0);
                                   if d2>d0 then
                                   begin
                                      d2:=d0;
                                      dn_max:=SpaceChargeArray[x,y,z].dn;
                                   end;
                                end;

                       for z := 0 to NodeCount_z-1  do
                         for x := 0 to NodeCount_x-1 do
                           for y := 0 to NodeCount_y-1 do
                              if SpaceChargeArray[x,y,z].Material=mSC then
                                begin
                                   d0:= SpaceChargeArray[x,y,z].dp;
                                   Newton_p(x,y,z);
                                   d0:=abs(SpaceChargeArray[x,y,z].dp-d0);
                                   if d3>d0 then
                                   begin
                                      d3:=d0;
                                      dp_max:=SpaceChargeArray[x,y,z].dp;
                                   end;
                                end;
                     end;
               until (m>0) and ((d1<=abs(dPhi_max*epsilon)) or (abs(dPhi_max)*PhiScale<1e-10))
                           and ((d2<=abs(dn_max*epsilon)) or (abs(dn_max)*NScale<1e-10))
                           and ((d3<=abs(dp_max*epsilon)) or (abs(dp_max)*NScale<1e-10));


              dw_k:=0;
              w_k:=0;

              for z := 0 to NodeCount_z-1  do
                for x := 0 to NodeCount_x-1 do
                  for y := 0 to NodeCount_y-1 do
                     begin
                          if (SpaceChargeArray[x,y,z].Material<>mMetal)and (not SpaceChargeArray[x,y,z].OhmicContact) then
                            begin
                              if SpaceChargeArray[x,y,z].Material=mSC then
                               begin
                                 SCIndex:=SpaceChargeArray[x,y,z].SCIndex;
                                 Offset:=Semiconductors[SCIndex].E_offset;
                               end
                              else
                                Offset:=0;
                              SpaceChargeArray[x,y,z].Phi:=SpaceChargeArray[x,y,z].Phi+SpaceChargeArray[x,y,z].dPhi;
                              dw_k:=dw_k+abs(SpaceChargeArray[x,y,z].dPhi)*PhiScale;
                              w_k:=w_k+abs(SpaceChargeArray[x,y,z].Phi+Offset)*PhiScale;
                            end;
                          if ((SpaceChargeArray[x,y,z].Material=mSC) and (not SpaceChargeArray[x,y,z].OhmicContact)) then
                            begin
                              if PreCalc then //Precalc (i.e. without Continuity equations)
                                begin
                                  SpaceChargeArray[x,y,z].n:=n(Semiconductors[SCIndex].E_f,SpaceChargeArray[x,y,z].Phi+Offset,Semiconductors[SCIndex].Inv_C,SCIndex);
                                  SpaceChargeArray[x,y,z].p:=p(Semiconductors[SCIndex].E_f,SpaceChargeArray[x,y,z].Phi+Offset,Semiconductors[SCIndex].Inv_V,SCIndex);
                                end
                              else
                                begin // (with continuity equations)
                                  SpaceChargeArray[x,y,z].n:=SpaceChargeArray[x,y,z].n+SpaceChargeArray[x,y,z].dn;
                                  SpaceChargeArray[x,y,z].p:=SpaceChargeArray[x,y,z].p+SpaceChargeArray[x,y,z].dp;
                                  //SpaceChargeArray[x,y,z].Phi_Gamma_N:=SpaceChargeArray[x,y,z].Phi_Gamma_N+0.1*omega*(Update_Phi_Gamma_n(x,y,z)-SpaceChargeArray[x,y,z].Phi_Gamma_N);
                                  //SpaceChargeArray[x,y,z].Phi_Gamma_P:=SpaceChargeArray[x,y,z].Phi_Gamma_P+0.1*omega*(Update_Phi_Gamma_p(x,y,z)-SpaceChargeArray[x,y,z].Phi_Gamma_P);
                                  SpaceChargeArray[x,y,z].Phi_Gamma_N:=Update_Phi_Gamma_n(x,y,z);
                                  SpaceChargeArray[x,y,z].Phi_Gamma_P:=Update_Phi_Gamma_p(x,y,z);

                                  dw_k:=dw_k+(abs(SpaceChargeArray[x,y,z].dn)+abs(SpaceChargeArray[x,y,z].dp))*Nscale;
                                  w_k:=w_k+(abs(SpaceChargeArray[x,y,z].n)+abs(SpaceChargeArray[x,y,z].p))*Nscale;
                                end;
                            end;
                     end;

               if PreCalc and (k>1) then
                 begin
                   PreCalc:=dw_k/w_k>dPhi_Grenz;//*1000;
                   if not PreCalc then
                     begin
                        Output('\n\nCalculation of system in thermal equilibrium done. Continuing with non-equilibrium.\n\n');
                        dw_k:=dw_k*10;
                        GetE_Fqn:=GetE_Fqn2;
                        GetE_Fqp:=GetE_Fqp2;
                        for x := 0 to NodeCount_x-1 do
                           for y := 0 to NodeCount_y-1 do
                               for z := 0 to NodeCount_z-1 do
                                  if (SpaceChargeArray[x,y,z].Material=mSC) and (not SpaceChargeArray[x,y,z].OhmicContact) then
                                    begin
                                       SpaceChargeArray[x,y,z].Phi_Gamma_N:=Update_Phi_Gamma_n(x,y,z);
                                       SpaceChargeArray[x,y,z].Phi_Gamma_P:=Update_Phi_Gamma_p(x,y,z);
                                    end;
                     end;
                 end;

               Phi_surf:=SpaceChargeArray[checkpoint_x,checkpoint_y,checkpoint_z].Phi;
               t_diff:=millisecondsbetween(t_start,now)*0.001;
               t_start:=now;
               if assigned(ProgressNotification) then
                 ProgressNotification(k,m,-Phi_surf,dw_k/w_k,t_diff);

           end;

      end;


end.

