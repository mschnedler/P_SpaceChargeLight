unit U_Common;

{$IFDEF FPC}
  {$MODE DELPHI}
{$ENDIF}

interface
uses SysUtils, Math;

const
      GoldenSection= 0.381966;
      Max_Int_Rekursion=200;
type
       TVec= record
         x,y,z:double;
       end;

       TStringDynArray = array of string;

       TFunctionParameters =  array of double;
       TCommonFunction     =  function(x: double; Parameters:TFunctionParameters): double;

       TMaterial = (mMetal, mSC, mVac);

       TBandGapType = (bgtDirect, bgtIndirect);

       TSpaceChargePoint= record
          Phi: double;
          Phi_Gamma_n:double;
          Phi_Gamma_p:double;
          Material: TMaterial;
          n: double;
          p: double;
          dPhi,dn,dp: double;
          dF1dPhi,dF2dn,dF3dp: double;
          Polarisation: TVec;          //Polarisation in 1/nm^2 (devided by e)
          SCIndex: byte;               //Semiconductor index
          SurfX,SurfY,SurfZ:boolean ;  //surface indicator
          IntX,IntY,IntZ:boolean;      //intface indicator
          SurfCount,IntCount: integer; //surface and interface Count
          divP: double;                //constant value of divergence P in 1/nm^3 (devided by e)
       end;

       TRhoList=record
         Value:     array of double;
         Count:     integer;
         Phi_Start: double;
         Phi_End:   double;
       end;

       // Properties of the semiconductor:
       TSemiconductor= record
          E_g:         double;       // Bandgap of semiconductor in eV
          E_offset:    double;       // Band edge offset between aligned semiconductors in eV
          N_D:         double;       // Donor concentration in 1/m³ (in Parameter.txt: 1/cm^3)
          N_A:         double;       // Acceptor concentration in 1/m³ (in Parameter.txt: 1/cm^3)
          N_C:         double;       // Effective carrier concentration of conduction band in 1/nm^3
          N_V:         double;       // Effective carrier concentration of valence band in 1/nm^3
          E_D:         double;       // Energy of donor level in eV
          E_A:         double;       // Energy of acceptor level in eV
          m_cb:        double;       // density of states effective mass of e- in conduction band as fraction of me
          m_vb:        double;       // density of states effective mass of h in valence band as fraction of me
          Epsi_s:      double;       // dielectric constant of semiconductor
          Chi:         double;       // electron affinity of semiconductor
          E_f:         double;       // Fermi level
          Interface_x: integer;      // right hand side interface of this semiconductor
          Const_Interface_charge: double;  // constant, fixed charge at the interface in 1/m^2
          Const_Surface_charge: double;    // constant, fixed charge at the surface in 1/m^2
          mu_n,mu_p:   double;       // Mobility of electrons and holes [m^2/(Vs)]
          n0,p0:       double;       // Electron and hole density (excluding light excited carriers in equilibrium) [1/nm^3]
          c_opt:       double;       // Electron / Hole density due to light excited carriers in equilibrium [1/nm^3]
          G_photo:     double;       // Generation rate of light excited carriers [1/(nm^3*s)]
          light_on:    boolean;      // Switch that "turns on the light".
          P_opt, E_ph,
          A_opt, alpha,
          tau:         double;       // Optical power of the laser [eV/s], photon energy [eV], Laser-spotsize [nm^2],
                                     // Absorption coefficient of semiconductor [1/nm], minority carrier lifetime of semiconductor [s]
          E_CNL:       double;       // Charge neutrality level  [eV]
          E_SS:        double;       // Peak position of the surface state distribution [eV]
                                     // NOTE:    E_CNL and E_SS are DEFINIED RELATIVE TO THE SURFACE VALENCE BAND EDGE,
                                     //          because the surface charge distribution moves with the
                                     //          surface band edges
          FWHM :       double;       // Width of the surface state distribution (full width half maximum) [eV]
          surface_charge_density:  double;    // Surface charge density of surface states [1/m^2]
          Dn,Dp:       double;       // Diffusion coefficient of electron and holes.
                                     // Common definition: [nm^2/s]
          C_Rate:      double;       // Constant for computational enhancement [nm^3/s]
          C_Poisson:   double;       // Constant for computational enhancement [V*nm]
          Epsi_semi:   double;       // Permittivity of the semiconductor in [C/(V*nm)]
          C_NA:        double;       // Acceptor concentration in [1/nm^3]
          C_ND:        double;       // Donor concentration in [1/nm^3]
          SurfRhoList: TRhoList;     // Used to store Surf_rho(Phi) values to speed up calculation
          n_List:      TRhoList;     // Used to store n(phi) values to speed up calculation
          p_List:      TRhoList;     // Used to store p(phi) values to speed up calculation
          Inv_V:       boolean;      // allow inversion with holes in valence band (n-type)
          Inv_C:       boolean;      // allow inversion with electrons in conduction band (p-type)
          BandGapType: TBandGapType; // Defines, whether we are dealing with an direct or indirect band gap

          Discontinuity_VB: double;  // Defines the valence band discontinuity relativ to Semiconductor 0.
          Discontinuity_CB: double;  // Defines the conduction band discontinuity relativ to Semiconductor 0.
          Polarisation: TVec;        // Defines the Polarisation field of the semiconductor [1/nm^2]
       end;

       TSpaceChargeArray= array of array of array of TSpaceChargePoint;

       function Integrate2(FunctionToIntegrate:TCommonFunction; Parameters:TFunctionParameters; From:double; Too:double;dx:double):double;
       function  Integrate(f:TCommonFunction; Parameters:TFunctionParameters; a:double; b:double; eps:double):double;
       function GoldenSectionInteration(f:TCommonFunction; Parameters:TFunctionParameters; var a:double; var b:double; eps:double):integer;
       procedure Output(s:string);
       procedure SplitString(Str:String; Delimiter: string; var DelStrAry: TStringDynArray);

var
       Semiconductors:      array of TSemiconductor;
       Semiconductor_Count: Integer;
       T:        double;       // Temperature in K
       d:        double;       // Tip-sample distance in m
       phi_m:    double;       // metal work function

       //Further variables used for calculation
       Tip_Radius,
       Tip_Apex_Angle:   double;   // in nm and degree
       C_kT:             double;   // Constant for computational enhancement 1/(kT) [1/eV]

       // Values for numerical evaluation:
        stepwidth_z:   double;   // width between two lateral points in the band bending iteration in m
                                // see also procedure i_tunnel and procedure phi
                                // Standard-value: 1e-10 m
       z_stop:        double;   // upper limit for z-integration in transmission coefficients.
                                // Standard-value: 200 nm
       I_V_Output:    boolean;  // is true, if the I(V) curve should be calculated and stored
       V_list:        array of double;  // Voltage list to be calculated
       I_V_filename:  string;   // filename for storage of I(V) curve
       BB_Z_output:   boolean;  // is true, if the band bending as a function of z should be calculated and stored for all voltage steps
       BB_Z_filename: string;   // filename for storage of the bend bending files (there should be a placeholder for the format command
                                // in the filename, like %.2e, see "format"-command in pascal help)
       BB_V_output:   boolean;  // is true, if the band bending as a function of v should be calculated and stored
       BB_V_filename: string;   // filename for storage of the bend bending file.
       EFQ_V_output:  boolean;  // is true, if the quasi-fermi levels and quasi-effective masses as a function of v should be calculated and stored
       EFQ_V_Filename: string;  // filename for storage of the quasi-fermilevels file
       Phi3D_output:  boolean;  // is true, if whole phi array should be saved
       Phi3D_filename:string;   // filename for storage of phi array. File extension will be added automatically

       Rho3D_output:     boolean;    // is true, if whole rho (n and p) array should be saved
       Rho3D_n_filename: string;     // filename for storage of n array. File extension will be added automatically
       Rho3D_p_filename: string;     // filename for storage of p array. File extension will be added automatically

       Rho_central_output: boolean;  // is true, if the n and p concentration as a function of z should be calculated and stored for all voltage steps
       Rho_central_filename: string; // filename for storage of the n and p concentration files (there should be a placeholder for the format command
                                     // in the filename, like %.2e, see "format"-command in pascal help)
       Interim_solution_Intervall:integer; //interval for storing an iterim solution (0= no storage)
       Interim_solution_x:integer;   // number of y-z plane, that should be stored
       Interim_solution_Path: string;// file path for storage of iterim solutions
       Interim_solution_n_filename: string; //filename for electron concentrations
       Interim_solution_p_filename: string; //filename for hole concentrations
       Interim_solution_phi_filename: string; //filename for phi array
       Interim_solution_ND_filename: string; //filename for ionized donor concentration
       Interim_solution_NA_filename: string; //filename for ionized acceptor concentration
       Interim_solution_Phi_Quantum_n_filename: string; //filename for quantum potential of electrons
       Interim_solution_Phi_Quantum_p_filename: string; //filename for quantum potential of holes
       E_FQ_C,E_FQ_V: double;        //Quasi Fermilevels for conduction band and valence band at the point opposite to the tip apex
       FjList:      TRhoList;
implementation


       procedure Output(s:string);
       begin
          s:=StringReplace(s,'\n',#13#10,[rfReplaceAll]);
          writeln(s);
       end;


       procedure SplitString(Str:String; Delimiter: string; var DelStrAry: TStringDynArray);
       var i:integer;
           n:integer;
       begin
         n:=0;
         setlength(DelStrAry,n);
         while pos(Delimiter,Str)> 0 do
           begin
             i:=pos(Delimiter,Str);
             inc(n);
             setlength(DelStrAry,n);
             DelStrAry[n-1]:=copy(Str,1,i-1);
             Str:=copy(Str,i+length(Delimiter),maxint);
           end;
         inc(n);
         setlength(DelStrAry,n);
         DelStrAry[n-1]:=Str;
       end;

       // *********************************************************************
       // Lyness's Modified Adaptive Simpson's method for integration
       // Journal of the ACM 16 483-495 (1969)
       // DOI: 10.1145/321526.321537
       // A numerical quadrature method that recursively bisects the interval
       // until the precision is high enough
       // See also:
       // rosettacode.org/wiki/Numerical_integration/Adaptive_Simpson%27s_method
       // **********************************************************************

       function Quad_Simpsons_mem(f:TCommonFunction; Parameters:TFunctionParameters; a,b,fa, fb: double; var m,fm:double): double;
       begin
         m:=(a+b)/2;
         fm:=f(m,parameters);
         result:=(b-a)/6*(fa+4*fm+fb);
       end;

       function Quad_Asr(f:TCommonFunction; Parameters:TFunctionParameters; recursion_level:integer; a, b, fa, fb, eps ,whole, m, fm: double): double;
       var lm,rm: double;
           left,right,flm,frm,r1,r2:double;
           delta:double;
       begin
         left:= Quad_Simpsons_mem(f,parameters, a, m, fa, fm,lm,flm);
         right:=Quad_Simpsons_mem(f,parameters, m, b, fm, fb,rm,frm);
         Delta:=left + right - whole;

         if ((abs(Delta)<=15* eps) and (recursion_level>2)) then
            begin
              result:=left + right + delta/15
            end
         else if recursion_level>Max_Int_Rekursion then
            begin
              result:=left + right + delta/15;
              //output('*** WARNING: Max. rekursion depth reached in Quad_Asr (adaptive simpson integration).');
            end
         else if (eps/2=eps) or (a=lm)then
            begin
              result:=left + right + delta/15;
              //output('*** WARNING: Interval refinement reached machines precision in Quad_Asr (adaptive simpson integration).');
            end
         else
           begin
            r1:= Quad_Asr(f,parameters,recursion_level+1,a,m,fa,fm,eps/2,left,lm,flm);
            r2:= Quad_Asr(f,parameters,recursion_level+1,m,b,fm,fb,eps/2,right,rm,frm);
            result:=r1+r2;
           end;
       end;

       function Integrate(f:TCommonFunction; Parameters:TFunctionParameters; a:double; b:double; eps:double):double;
       var fa,fb:Double;
           res,fm,whole: double;
           m:double;
       begin
         if (a=b) then
          begin
            result:=0;
            exit;
          end;
         fa:=f(a,parameters);
         fb:=f(b,parameters);
         whole:= Quad_Simpsons_mem(f,parameters,a,b,fa,fb,m,fm);
         result:=Quad_Asr(f,parameters,0,a,b,fa,fb,eps,whole,m,fm);
       end;

       function Integrate2(FunctionToIntegrate:TCommonFunction; Parameters:TFunctionParameters; From:double; Too:double;dx:double):double;
       var width,dummy:double;
           i,steps:cardinal;
           vz:integer;
       begin
          dx:=0.0001;
          if (dx=0) or (Too=From) then
            begin
              result:=0;
              exit;
            end;
          width:=Too-From;
          dx:=abs(dx)*(Too-From)/abs(Too-From);
          if width/dx>high(steps) then
             begin
               output('Integration-Error: dx is set too small. Using smallest acceptable dx = '+floattostr(width/high(steps)));
               dx:=width/high(steps);
             end;
          steps:=round(width/dx);
          result:=0;

          for i := 1 to steps do
              result:=result+FunctionToIntegrate(From+(i-0.5)*dx, Parameters);

          result:=result*dx;
       end;





       function GoldenSectionInteration(f:TCommonFunction; Parameters:TFunctionParameters; var a:double; var b:double; eps:double):integer;
       var Tau:double;
           range,range_old,t:double;
           x1,x2:double;
           fx1,fx2:Double;
           Num_Iter:integer;
       begin
         result:=0;
         if a=b then exit;
         if a>b then
           begin
             t:=a;
             a:=b;
             b:=t;
           end;
         range:=b-a;
         x1:=a+GoldenSection*range;
         x2:=b-GoldenSection*range;
         fx1:=f(x1,parameters);
         fx2:=f(x2,parameters);
         range_old:=0;
         Num_Iter:=0;
         while (range>eps) and (range_old<>range) do
          begin
            inc(Num_Iter);
            range_old:=range;
            if fx1<fx2 then
              begin
                  b:=x2;
                  x2:=x1;
                  fx2:=fx1;
                  range:=b-a;
                  x1:=a+GoldenSection*range;
                  fx1:=f(x1,parameters);
               end
              else
               begin
                  a:=x1;
                  range:=b-a;
                  x1:=x2;
                  fx1:=fx2;
                  x2:=b-GoldenSection*range;
                  fx2:=f(x2,Parameters);
               end;
          end;
          result:=Num_Iter;
          if (abs(b-a)>eps) then
               raise EMathError.Create(format('Golden section iteration did not converge:\n'+
                                       'E_F not found within %d steps and precision = %.5e.\n'+
                                       'Try to increase decrease precision!',
                                       [Num_Iter,eps]));
       end;



end.
