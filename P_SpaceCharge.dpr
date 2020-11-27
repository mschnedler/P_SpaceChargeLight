program P_SpaceCharge;

{$IFDEF FPC}
  {$MODE DELPHI}
{$ELSE}
  {$APPTYPE CONSOLE}
  {$R *.res}
{$ENDIF}


uses
  SysUtils,
  Math,
  Classes,
  DateUtils,
  U_Constants in 'U_Constants.pas',
  U_Common in 'U_Common.pas',
  U_CarrierConcentration in 'U_CarrierConcentration.pas',
  U_TunnelCurrent_2 in 'U_TunnelCurrent_2.pas',
  U_FDMSolver in 'U_FDMSolver.pas';

type
       TGridPointArray=array of double;
var
       // Definition of variables:
       PhiScale_tip:                          double;

       Voltage:                               double;   // current tip-sample voltage step
       Parameterlist: array of string;                  // Paremeters read from parameter input file
       x0,y0,z0:                              integer;  // Parameters for Grid-Transformation-Function
       dx_min,dy_min,dz_min:                  double;   // Parameter for Grid-Transformation-Function
       NodeCount_x_Start, NodeCount_y_Start,            // Number of Points in x,y,z direction of the space charge grid
       NodeCount_z_Start :                    integer;  // before the first refinement is performed

       refinement_steps:                      integer;   // Defines the number of grid-node refinements.
                                                         // Each refinement doubles the number of grid points in each direction.
       Ignore_1D_estimation:                  boolean;   // Switch defining, whether there should be a warning, if the
                                                         // lateral extension of 1D estimation of the band bending is larger
                                                         // the size of the grid.
       Quiet:                            boolean;        // Flag indicating, if user inputs are allowed
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

       // Read parameter file and set variables
       function ReadParameterFile(filename:string):boolean;
       var tf                 : textfile;
            linearray          : array of string;
            numlines           : integer;
            numparameters      : integer;
            assigned_parameters: array of boolean;
            i,j,k              : integer;
            paramname          : string;
            value              : string;
            line               : string;
            sda                : TStringDynArray;
       begin
         result:=false;
         if not fileexists(filename) then
           begin
            output('file not found: '+filename);
            exit;
           end;
         numparameters := length(ParameterNames);
         setlength(assigned_parameters, numparameters);
         setlength(Parameterlist, numparameters);
         for i := 0 to numparameters-1 do
            assigned_parameters[i]:=false;
         assignfile(tf,filename);
         reset(tf);
         numlines:=0;
         while not eof(tf) do
           begin
             inc(numlines);
             setlength(linearray,numlines);
             readln(tf,linearray[numlines-1]);
           end;
         closefile(tf);
         output('\nRead '+inttostr(numlines)+' lines from parameter file '+filename+':');
         for i := 0 to numlines-1 do
           begin
              if pos('#',linearray[i])=0 then
                line:=linearray[i]
              else
                line:=copy(linearray[i],1,pos('#',linearray[i])-1);
              if pos('=',line)>0 then
               begin
                  paramname := trim( copy(line,1,pos('=',line)-1 ));
                  value     := trim( copy(line,pos('=',line) + 1,maxint));
                  for j:= 0 to numparameters-1 do
                    begin
                      if uppercase(ParameterNames[j])=uppercase(paramname) then
                        begin
                          ParameterList[j]:=value;
                          assigned_parameters[j]:=true;
                          output(ParameterNames[j]+':'+#9+ value);
                        end;
                    end;
               end;
           end;
          result:=true;
          for i := 0 to numparameters-1 do
            begin
              result:=result and assigned_parameters[i];
              if not assigned_parameters[i] then
                output('missing parameter in input file: '+ParameterNames[i]);
            end;
          i:=0;
          if result then
          begin
             try
              T      := strtofloat(ParameterList[i]);
              inc(i);
              d      := strtofloat(ParameterList[i]);
              inc(i);
              phi_m  := strtofloat(ParameterList[i]);
              inc(i);
              // Set parameters of semiconductor:
              SplitString(ParameterList[i],',',sda);
              Semiconductor_Count:=length(sda);
              setlength(Semiconductors,Semiconductor_Count);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Epsi_s:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Chi:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].m_cb:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].m_vb:=strtofloat(sda[k]);
              inc(i);
              //m_vb=exp(2.*alog(sqrt(heavyhole_effectivemass**3)+sqrt(lighthole_effectivemass**3))/3.)
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].E_g:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].N_D:=strtofloat(sda[k]) * power(100,3);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].N_A:=strtofloat(sda[k]) * power(100,3);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].E_D:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].E_A:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Inv_V:=strtobool(trim(sda[k]));
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Inv_C:=strtobool(trim(sda[k]));
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Polarisation.x:=strtofloat(trim(sda[k]));
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Polarisation.y:=strtofloat(trim(sda[k]));
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Polarisation.z:=strtofloat(trim(sda[k]));
              inc(i);
              // Set parameters of numerical evaluation:
              stepwidth_z    := strtofloat(ParameterList[i]);
              inc(i);
              I_V_output     := strtobool(ParameterList[i]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              setlength(V_list,length(sda));
              for k := 0  to length(sda)-1 do
                V_list[k]:=strtofloat(sda[k]);
              inc(i);
              I_V_filename   := ParameterList[i];
              inc(i);
              BB_Z_output    := strtobool(ParameterList[i]);
              inc(i);
              BB_Z_filename  := ParameterList[i];
              inc(i);
              BB_V_output    := strtobool(ParameterList[i]);
              inc(i);
              BB_V_filename  := ParameterList[i];
              inc(i);
              Phi3D_output  := strtobool(ParameterList[i]);
              inc(i);
              Phi3D_filename:= ParameterList[i];
              inc(i);
              NodeCount_x_Start:=strtoint(ParameterList[i]);
              inc(i);
              NodeCount_y_Start:=strtoint(ParameterList[i]);
              inc(i);
              NodeCount_z_Start:=strtoint(ParameterList[i]);
              inc(i);
              x_length:=strtofloat(ParameterList[i]);
              inc(i);
              y_length:=strtofloat(ParameterList[i]);
              inc(i);
              z_length:=strtofloat(ParameterList[i]);
              inc(i);
              dx_min:=strtofloat(ParameterList[i]);
              inc(i);
              dy_min:=strtofloat(ParameterList[i]);
              inc(i);
              dz_min:=strtofloat(ParameterList[i]);
              inc(i);
              dPhi_Grenz:=strtofloat(ParameterList[i]);
              inc(i);
              omega:=strtofloat(ParameterList[i]);
              inc(i);
              Tip_Radius:=strtofloat(ParameterList[i]);
              inc(i);
              Tip_Apex_Angle:=strtofloat(ParameterList[i]);
              inc(i);
              refinement_steps:=strtoint(ParameterList[i]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].E_SS:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].E_CNL:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].FWHM:=strtofloat(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].surface_charge_density:=strtofloat(sda[k])*power(100,2); // 1/m^2;
              inc(i);
              z_stop:=strtofloat(ParameterList[i]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].mu_n:=strtofloat(sda[k])/1e4; //Parameterfile: cm^2/Vs, hier: m^2/Vs
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].mu_p:=strtofloat(sda[k])/1e4; //Parameterfile: cm^2/Vs, hier: m^2/Vs
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].P_opt:=strtofloat(sda[k])/eVToJ;   //Parameterfile: W, hier: eV/s
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].E_ph:=strtofloat(sda[k]); //eV
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].A_opt:=strtofloat(sda[k]); //nm^2
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].alpha:=strtofloat(sda[k])*1e-7; //Parameterfile: 1/cm, hier: 1/nm
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].tau:=strtofloat(sda[k]); //s
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].light_on:=strtobool(trim(sda[k]));
              inc(i);
              Rho3D_output:=strtobool(ParameterList[i]);
              inc(i);
              Rho3D_n_filename:=ParameterList[i];
              inc(i);
              Rho3D_p_filename:=ParameterList[i];
              inc(i);
              Rho_central_output:=strtobool(ParameterList[i]);
              inc(i);
              Rho_central_filename:=ParameterList[i];

              inc(i);
              Interim_solution_Intervall:=strtoint(ParameterList[i]);
              inc(i);
              Interim_solution_x:=strtoint(ParameterList[i]);
              inc(i);
              Interim_solution_Path:=ParameterList[i];
              inc(i);
              Interim_solution_n_filename:=ParameterList[i];
              inc(i);
              Interim_solution_p_filename:=ParameterList[i];
              inc(i);
              Interim_solution_phi_filename:=ParameterList[i];
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Interface_x:=strtoint(sda[k]);
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                if uppercase(trim(sda[k]))='DIRECT' then
                   Semiconductors[k].BandGapType:=bgtDirect
                else if uppercase(trim(sda[k]))='INDIRECT' then
                   Semiconductors[k].BandGapType:=bgtIndirect
                else
                   raise Exception.Create(ErrorMessages[4]);
              inc(i);
              Interim_solution_ND_filename:=ParameterList[i];
              inc(i);
              Interim_solution_NA_filename:=ParameterList[i];
              inc(i);
              Interim_solution_phi_quantum_n_filename:=ParameterList[i];
              inc(i);
              Interim_solution_phi_quantum_p_filename:=ParameterList[i];
              inc(i);
              EFQ_V_output:=strtobool(ParameterList[i]);
              inc(i);
              EFQ_V_filename:=ParameterList[i];
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Const_Interface_charge:=strtofloat(sda[k])*100*100; //1/m^2
              inc(i);
              SplitString(ParameterList[i],',',sda);
              if Semiconductor_Count<>length(sda) then
                raise Exception.Create(ErrorMessages[0]);
              for k := 0  to Semiconductor_Count-1 do
                Semiconductors[k].Const_surface_charge:=strtofloat(sda[k])*100*100; //1/m^2
             except
                  on E : Exception do
                    begin
                     output('\nErroneous parameter: '+ParameterNames[i]+' -> Value "'+ParameterList[i]+'" is not allowed.');
                     output('Message: '+E.Message+'\n');
                     result:=false;
                    end;
             end;
          end;
       end;

      // Grid transformation function (from grid-indices to physical values)
      // Calculates the grid spacing such, that a maximal number of grid points
      // around 0 exhibit the minimal point separation dx (given by the user
      // in the parameter file as dx_min, dy_min, dz_min). Further, by doubling the
      // grid separation for each point it is ensured, that the grid points cover the
      // total length x_len of the simulation region.
      // Note: The user-input x_length, y_length, z_length is NOT the precise length
      //       in each of the three directions. It is only a lower limit. However
      //       due to the doubling of the grid spacing, the real size may differ from
      //       the user input by a factor of nearly 2.
      // Note: The numerical error of the calculation depends on the separation of
      //       the grid points. Larger grid spacing means larger numerical errors.
      //       See [2], p. 159 for error approximation.

      function GridToPhysical4(NodeCount:integer;dx,x_len:double;x,x0:integer):double;
      var i:integer;
          dxx:double;
          grenz:integer;
      begin
         result:=0;
         dxx:=dx;
         grenz:=NodeCount div 2;
         repeat
           result:=0;
           dxx:=dx;
           grenz:=grenz-1;
           for i := 0 to abs(NodeCount div 2)-1 do
             begin
                if (i>grenz) and ((i-grenz) mod 3 =0) then
                  dxx:=2*dxx;
                result:=result+dxx;
             end;
         until abs(result)>x_len/2;

         result:=0;
         dxx:=dx;

         for i := 0 to abs(x-x0)-1 do
           begin
              if (i>grenz) and ((i-grenz) mod 3 =0) then
                dxx:=2*dxx;
              result:=result+dxx*sign(x-x0);
           end;
      end;

      // Conversion of 3D grid indices x,y,z to physical coordinates px,py,pz
      procedure GridToPhysical_Array(x,y,z:integer; var px,py,pz:double);
      begin
        if (x=-1) or (x=NodeCount_x) then px:=0
          else
        px:=px_array[x];

        if (y=-1) or (y=NodeCount_y) then py:=0
          else
        py:=py_array[y];

        if (z=-1) or (z=NodeCount_z) then pz:=0
          else
        pz:=pz_array[z];
      end;

      // Conversion of a physical coordinate p to an integer grid index.
      function PhysicalToGrid(var ar:array of double; p:double):integer;
      var diff_array: array of double;
          i: Integer;
          amin:integer;
      begin
        setlength(diff_array,length(ar));
        for i := 0 to length(ar)-1 do
           diff_array[i]:=abs(ar[i]-p);
        amin:=0;

        for i := 1 to length(ar)-1 do
           if diff_array[i]<diff_array[amin] then amin:=i;

        result:=amin;
      end;

      // Finding the physical grid values g1,g2 (and the corresponding
      // grid indices i1,i2) that are located closest to an arbitrary physical
      // value p.
      function PhysicalToDiscretePhysical(var ar:array of double; p:double;var g1,g2:double; var i1,i2:integer):boolean;
      var i: Integer;
      begin
        i:=0;
        g2:=ar[0];
        g1:=ar[0]+(ar[1]-ar[0])*(-1);
        result:=true;
        while g2<p do
          begin
            inc(i);
            if i>=length(ar) then
               begin
                 g2:=ar[length(ar)-1]+(ar[length(ar)-1]-ar[length(ar)-2])*(-1) ;
                 g1:=ar[length(ar)-1];
                 i2:=length(ar);
                 i1:=length(ar)-1;
                 result:=false;
                 exit;
               end;
            g2:=ar[i];
            g1:=ar[i-1];
          end;
        if g2=p then g1:=g2;
        i2:=i;
        if g2=p then i1:=i2 else i1:=i2-1;
      end;

      // Initializing the arrays px_array, py_array, and pz_array that store
      // the physical positions of the mesh points. Additionally, the
      // point in z-direction, that is located closest to the foremost
      // tip-apex is set to the precise value of the tip-sample separation.
      // With the help of this method, the tip-sample separation can be
      // adjusted with higher precision as compared to the minimal grid-spacing.
      function Init_Coordinate_System:boolean;
      var x,y,z:integer;
      begin
        result:=true;
        if ceil(log2(x_length/2/dx_min))>NodeCount_x div 2 then
        begin
          output('Error: dx_min is set too small/ x_length is set too large.');
          result:=false;
          exit;
        end;
        if ceil(log2(y_length/2/dy_min))>NodeCount_y div 2 then
        begin
          output('Error: dy_min is set too small/ y_length is set too large.');
          result:=false;
          exit;
        end;
        if ceil(log2(z_length/2/dz_min))>NodeCount_z div 2 then
        begin
          output('Error: dz_min is set too small/ z_length is set too large.');
          result:=false;
          exit;
        end;
        setlength(px_array,NodeCount_x);
        setlength(py_array,NodeCount_y);
        setlength(pz_array,NodeCount_z);
        for x := 0 to NodeCount_x-1 do
          px_array[x]:=GridToPhysical4(NodeCount_x,dx_min,x_length,x,x0);// GridToPhysical(ax,bx,x0,x);
        for y := 0 to NodeCount_y-1 do
          py_array[y]:=GridToPhysical4(NodeCount_y,dy_min,y_length,y,y0); //GridToPhysical(ay,by,y0,y);
        for z := 0 to NodeCount_z-1 do
          pz_array[z]:=GridToPhysical4(NodeCount_z,dz_min,z_length,z,z0); //GridToPhysical(az,bz,z0,z);

        pz_array[PhysicalToGrid(pz_array,-d)]:=-d;
      end;

      // 3D Hyperbolic function in order to derive the shape of a
      // hyperbolic tip. Any other tip shape can be defined within this procedure
      // r:      tip-radius
      // alpha:  tip apex angle
      // d:      tip-sample distance
      function Tip(x,y:double; r, alpha, d:double):double;
      var e,a,b:double;
      begin
        b:=r/tan(alpha/360*2*Pi);
        a:=b*b/r;
        e:=sqrt(a*a+b*b);
        result:=a*sqrt( 1+ (x*x+y*y)/(b*b));
        result:=result-a+d;
      end;

    // Calculates a one-dimensional estimation of the screening length of the
    // electrostatic potential as a function of the potential at the surface of
    // the semiconductor and of the dopand concentrations.
    // In principle, this is the Debye length, where the thermal voltage is
    // replaced by the applied electrostatic potential.
    // This estimation is in analogy to [3]
    function Calc_1D_Estimation(aPhi:double;SC:integer):double;
     begin
       //Schottky barrier width
       //aPhi: Potential at the surface of the semiconductor
       if (Semiconductors[SC].N_D-Semiconductors[SC].N_A)<>0 then
         result:=sqrt(2*Epsi_0*Semiconductors[SC].Epsi_s*abs(aPhi)/
                      (abs(Semiconductors[SC].N_D-Semiconductors[SC].N_A)*e)
                     )*1e9
       else
         result:=1e10;
     end;

     function GetSemiconductorAtPosition(x:integer):integer;
     var i:integer;
     begin
         for i := 0 to Semiconductor_Count-1 do
           if Semiconductors[i].Interface_x > x then break;
         result:=max(i-1,0);
     end;

     // Find the semiconductor that is at the tip-apex
     // Note the vacuum energy of the semiconductor beneath the tip
     // is used as zero potential energy in the finite difference
     // scheme.
     function GetSemiconductorBelowTipApex:integer;
     begin
         result:=GetSemiconductorAtPosition(NodeCount_x_Start div 2);
     end;


      // Initialization of some constants and the so called
      // three-dimensional space charge array, that contains the
      // electrostatic potential, the carrier concentrations for holes and
      // electrons, and the type of material (Semiconductor, Vacuum, Tip)
      // for each point of the computation grid.
      procedure Init_SpaceChargeArray(Voltage:double);
      var x,y,z,SC: Integer;
          z_start: Integer;
          px,py,pz:double;
          Phi_tip: double;
          pz_start_min:double;
      begin
         pz_start_min:=10000;
         Phi_tip:=-(Voltage+Semiconductors[0].E_F-Semiconductors[0].E_g+(Phi_m-Semiconductors[0].Chi));

         output(format('\nPhi_tip = %.5e V',[-Phi_tip]));

         PhiScale_tip:=1/abs(Phi_tip);
         for y := 0 to NodeCount_y-1 do
           for x := 0 to NodeCount_x-1 do
              begin
                 GridToPhysical_Array(x,y,0,px,py,pz);
                 pz:=-Tip(px,py,Tip_Radius,Tip_Apex_Angle,d);

                 if abs(pz)>=z_length/2 then
                    z_start:=-1
                  else
                    z_start:=PhysicalToGrid(pz_array,pz);

                 if z_start>-1 then
                    pz_start_min:=min(-pz_array[z_start],pz_start_min);

                 for z :=  z_start downto 0 do
                        AssignMetal(x,y,z,Phi_tip);

                 for z := max( z_start+1,0) to NodeCount_z-1 do
                   begin
                     if z>=NodeCount_z div 2 then
                       begin
                         SC:=GetSemiconductorAtPosition(x);
                         AssignSC(x,y,z,SC)
                       end
                     else
                       AssignVacuum(x,y,z);
                   end;
              end;

         output('z_start_min = '+floattostr(pz_start_min)+' nm ('+inttostr(PhysicalToGrid(pz_array,pz_start_min)-NodeCount_z div 2-1)+' node points)');
         Checkpoint_x:=NodeCount_x div 2;
         Checkpoint_y:=NodeCount_y div 2;
         Checkpoint_z:=NodeCount_z div 2;
      end;



      // Initialization of Coordinate System, Space Charge Array (containing the
      // electrostatic potential phi and the charge carrier densities n,p)
      // Calculation of the 1D estimation of the electrostatic potential.
      // Resetting of the values NodeCount_x,NodeCount_y, and NodeCount_z to their
      // initial values in case, that one uses a Grid Refinement, that increases
      // these values during the calculation.
      function Init(Voltage:double):boolean;
      var  a:string;
           i:integer;
      begin


            SetNodeCount(NodeCount_x_Start,NodeCount_y_Start,NodeCount_z_Start);

            for i:=0 to Semiconductor_Count-1 do
              output(format('\n1D estimation of depletion region of semiconductor '+inttostr(i)+': %.5f nm\n',[Calc_1D_Estimation(max(1,abs(Voltage)),i)]));

            for i:=0 to Semiconductor_Count-1 do
              if (not Ignore_1D_estimation) and (not Quiet) and (Calc_1D_Estimation(max(1,abs(Voltage)),i)>z_length/2) then
              begin
                output('*** WARNING: physical z-length in computation is smaller than the 1D estimation.\nThis may result in a false solution of the poisson equation.\n\nContinue anyway [(y)es| (n)o| (a)lways yes]?');
                readln(a);
                if uppercase(a)='A' then
                  Ignore_1D_estimation:=true
                else if uppercase(a)<>'Y' then
                  begin
                    result:=false;
                    exit;
                  end;
                break;
              end;

            x0:=NodeCount_x div 2;
            y0:=NodeCount_y div 2;
            z0:=NodeCount_z div 2;

            result:=Init_Coordinate_System;
            if not result then exit;

            Init_SpaceChargeArray(Voltage);
      end;

      //Helping function for tri-linear interpolation of phi, n, and p
      function rF(i:integer;x1,x2,xn:double):double;
      begin
         if i=1 then
           if x1=xn then result:=1
           else if x2=xn then result:=0
                else result:=(x2-xn)/(x2-x1);

         if i=2 then
           if x1=xn then result:=0
           else if x2=xn then result:=1
                else result:=(xn-x1)/(x2-x1);
      end;

      // Interpolation of Phi on physical positions between node-points
      // This procedure is used in the refinement of the grid.
      function Interpolate_Trilinear_Phi(px,py,pz:double):double;
      var ix1,iy1,iz1,ix2,iy2,iz2:integer;
          x1,x2,y1,y2,z1,z2:double;
      begin
               PhysicalToDiscretePhysical(px_array,px,x1,x2,ix1,ix2);
               PhysicalToDiscretePhysical(py_array,py,y1,y2,iy1,iy2);
               PhysicalToDiscretePhysical(pz_array,pz,z1,z2,iz1,iz2);
               if ix1<0 then ix1:=0;
               if iy1<0 then iy1:=0;
               if iz1<0 then iz1:=0;
               if ix1>NodeCount_x-1 then ix1:=NodeCount_x-1;
               if iy1>NodeCount_y-1 then iy1:=NodeCount_y-1;
               if iz1>NodeCount_z-1 then iz1:=NodeCount_z-1;
               if ix2<0 then ix2:=0;
               if iy2<0 then iy2:=0;
               if iz2<0 then iz2:=0;
               if ix2>NodeCount_x-1 then ix2:=NodeCount_x-1;
               if iy2>NodeCount_y-1 then iy2:=NodeCount_y-1;
               if iz2>NodeCount_z-1 then iz2:=NodeCount_z-1;
               result:=         Phi(ix1,iy1,iz1)*rF(1,x1,x2,px)*rF(1,y1,y2,py)*rF(1,z1,z2,pz)
                               +Phi(ix2,iy1,iz1)*rF(2,x1,x2,px)*rF(1,y1,y2,py)*rF(1,z1,z2,pz)
                               +Phi(ix1,iy2,iz1)*rF(1,x1,x2,px)*rF(2,y1,y2,py)*rF(1,z1,z2,pz)
                               +Phi(ix2,iy2,iz1)*rF(2,x1,x2,px)*rF(2,y1,y2,py)*rF(1,z1,z2,pz)
                               +Phi(ix1,iy1,iz2)*rF(1,x1,x2,px)*rF(1,y1,y2,py)*rF(2,z1,z2,pz)
                               +Phi(ix2,iy1,iz2)*rF(2,x1,x2,px)*rF(1,y1,y2,py)*rF(2,z1,z2,pz)
                               +Phi(ix1,iy2,iz2)*rF(1,x1,x2,px)*rF(2,y1,y2,py)*rF(2,z1,z2,pz)
                               +Phi(ix2,iy2,iz2)*rF(2,x1,x2,px)*rF(2,y1,y2,py)*rF(2,z1,z2,pz);
      end;

      // Interpolation of n on physical positions between node-points
      // This procedure is used in the refinement of the grid.
      function Interpolate_Trilinear_n(px,py,pz:double):double;
      var ix1,iy1,iz1,ix2,iy2,iz2:integer;
          x1,x2,y1,y2,z1,z2:double;
      begin
               PhysicalToDiscretePhysical(px_array,px,x1,x2,ix1,ix2);
               PhysicalToDiscretePhysical(py_array,py,y1,y2,iy1,iy2);
               PhysicalToDiscretePhysical(pz_array,pz,z1,z2,iz1,iz2);
               if ix1<0 then ix1:=0;
               if iy1<0 then iy1:=0;
               if iz1<0 then iz1:=0;
               if ix1>NodeCount_x-1 then ix1:=NodeCount_x-1;
               if iy1>NodeCount_y-1 then iy1:=NodeCount_y-1;
               if iz1>NodeCount_z-1 then iz1:=NodeCount_z-1;
               if ix2<0 then ix2:=0;
               if iy2<0 then iy2:=0;
               if iz2<0 then iz2:=0;
               if ix2>NodeCount_x-1 then ix2:=NodeCount_x-1;
               if iy2>NodeCount_y-1 then iy2:=NodeCount_y-1;
               if iz2>NodeCount_z-1 then iz2:=NodeCount_z-1;
               result:=         Get_n(ix1,iy1,iz1)*rF(1,x1,x2,px)*rF(1,y1,y2,py)*rF(1,z1,z2,pz)
                               +Get_n(ix2,iy1,iz1)*rF(2,x1,x2,px)*rF(1,y1,y2,py)*rF(1,z1,z2,pz)
                               +Get_n(ix1,iy2,iz1)*rF(1,x1,x2,px)*rF(2,y1,y2,py)*rF(1,z1,z2,pz)
                               +Get_n(ix2,iy2,iz1)*rF(2,x1,x2,px)*rF(2,y1,y2,py)*rF(1,z1,z2,pz)
                               +Get_n(ix1,iy1,iz2)*rF(1,x1,x2,px)*rF(1,y1,y2,py)*rF(2,z1,z2,pz)
                               +Get_n(ix2,iy1,iz2)*rF(2,x1,x2,px)*rF(1,y1,y2,py)*rF(2,z1,z2,pz)
                               +Get_n(ix1,iy2,iz2)*rF(1,x1,x2,px)*rF(2,y1,y2,py)*rF(2,z1,z2,pz)
                               +Get_n(ix2,iy2,iz2)*rF(2,x1,x2,px)*rF(2,y1,y2,py)*rF(2,z1,z2,pz);
      end;

      // Interpolation of p on physical positions between node-points
      // This procedure is used in the refinement of the grid.
      function Interpolate_Trilinear_p(px,py,pz:double):double;
      var ix1,iy1,iz1,ix2,iy2,iz2:integer;
          x1,x2,y1,y2,z1,z2:double;
      begin
               PhysicalToDiscretePhysical(px_array,px,x1,x2,ix1,ix2);
               PhysicalToDiscretePhysical(py_array,py,y1,y2,iy1,iy2);
               PhysicalToDiscretePhysical(pz_array,pz,z1,z2,iz1,iz2);
               if ix1<0 then ix1:=0;
               if iy1<0 then iy1:=0;
               if iz1<0 then iz1:=0;
               if ix1>NodeCount_x-1 then ix1:=NodeCount_x-1;
               if iy1>NodeCount_y-1 then iy1:=NodeCount_y-1;
               if iz1>NodeCount_z-1 then iz1:=NodeCount_z-1;
               if ix2<0 then ix2:=0;
               if iy2<0 then iy2:=0;
               if iz2<0 then iz2:=0;
               if ix2>NodeCount_x-1 then ix2:=NodeCount_x-1;
               if iy2>NodeCount_y-1 then iy2:=NodeCount_y-1;
               if iz2>NodeCount_z-1 then iz2:=NodeCount_z-1;
               result:=         Get_p(ix1,iy1,iz1)*rF(1,x1,x2,px)*rF(1,y1,y2,py)*rF(1,z1,z2,pz)
                               +Get_p(ix2,iy1,iz1)*rF(2,x1,x2,px)*rF(1,y1,y2,py)*rF(1,z1,z2,pz)
                               +Get_p(ix1,iy2,iz1)*rF(1,x1,x2,px)*rF(2,y1,y2,py)*rF(1,z1,z2,pz)
                               +Get_p(ix2,iy2,iz1)*rF(2,x1,x2,px)*rF(2,y1,y2,py)*rF(1,z1,z2,pz)
                               +Get_p(ix1,iy1,iz2)*rF(1,x1,x2,px)*rF(1,y1,y2,py)*rF(2,z1,z2,pz)
                               +Get_p(ix2,iy1,iz2)*rF(2,x1,x2,px)*rF(1,y1,y2,py)*rF(2,z1,z2,pz)
                               +Get_p(ix1,iy2,iz2)*rF(1,x1,x2,px)*rF(2,y1,y2,py)*rF(2,z1,z2,pz)
                               +Get_p(ix2,iy2,iz2)*rF(2,x1,x2,px)*rF(2,y1,y2,py)*rF(2,z1,z2,pz);
      end;

      // Grid refinement: If a solution for a small number of grid points was found,
      //                  the values of Phi, n, and p can be interpolated and stored
      //                  into a grid exhibiting a significant larger number of
      //                  grid points. This refines grid is used as initial solution
      //                  for a second run of the finite difference iteration scheme.
      // Note:            Slight differences between the surface values of phi, n, p
      //                  usually arise from (a) a changed grid spacing and hence a
      //                  changed numerical error of the calculation and (b) due to
      //                  a slightly different tip-shape, that is approximated with
      //                  with higher resolution, if the number of grid points is
      //                  increased.
      procedure RefineGrid(New_NodeCount_x, New_NodeCount_y, New_NodeCount_z:integer);
      var tmp: TSpaceChargeArray;
          x0n,y0n,z0n:integer;
          px_array_n,py_array_n,pz_array_n: Array of double;
          dx_min_n,dy_min_n,dz_min_n:double;
          xn,yn,zn:double;
          i,j,k: Integer;
      begin
         output('Refining grid to: '+inttostr(New_NodeCount_x)+'x'+inttostr(New_NodeCount_y)+'x'+inttostr(New_NodeCount_z)+' node points.');
         setlength(tmp,New_NodeCount_x,New_NodeCount_y,New_NodeCount_z);

         dx_min_n:=dx_min;///2;
         dy_min_n:=dy_min;///2;
         dz_min_n:=dz_min;///2;

         x0n:=New_NodeCount_x div 2;
         y0n:=New_NodeCount_y div 2;
         z0n:=New_NodeCount_z div 2;

         setlength(px_array_n,New_NodeCount_x);
         setlength(py_array_n,New_NodeCount_y);
         setlength(pz_array_n,New_NodeCount_z);

         for i := 0 to New_NodeCount_x-1 do
           px_array_n[i]:=GridToPhysical4(New_NodeCount_x,dx_min_n,x_length,i,x0n); //(axn,bxn,x0n,i);

         for i := 0 to New_NodeCount_y-1 do
           py_array_n[i]:=GridToPhysical4(New_NodeCount_y,dy_min_n,y_length,i,y0n);//GridToPhysical2(ayn,byn,y0n,i);

         for i := 0 to New_NodeCount_z-1 do
            pz_array_n[i]:=GridToPhysical4(New_NodeCount_z,dz_min_n,z_length,i,z0n);//GridToPhysical2(azn,bzn,z0n,i);

         pz_array_n[PhysicalToGrid(pz_array_n,-d)]:=-d;

         for k := 0 to New_NodeCount_z-1 do
           for j := 0 to New_NodeCount_y-1 do
              for i := 0 to New_NodeCount_x-1 do
                   begin
                     xn:=GridToPhysical4(New_NodeCount_x,dx_min_n,x_length,i,x0n);
                     yn:=GridToPhysical4(New_NodeCount_y,dy_min_n,y_length,j,y0n);
                     zn:=GridToPhysical4(New_NodeCount_z,dz_min_n,z_length,k,z0n);
                     tmp[i,j,k].Phi:=Interpolate_Trilinear_Phi(xn,yn,zn);
                     tmp[i,j,k].n:=Interpolate_Trilinear_n(xn,yn,zn);
                     tmp[i,j,k].p:=Interpolate_Trilinear_p(xn,yn,zn);
                   end;

         setlength(px_array,New_NodeCount_x);
         setlength(py_array,New_NodeCount_y);
         setlength(pz_array,New_NodeCount_z);

         for i := 0 to New_NodeCount_x-1 do
           px_array[i]:=px_array_n[i];

         for i := 0 to New_NodeCount_y-1 do
           py_array[i]:=py_array_n[i];

         for i := 0 to New_NodeCount_z-1 do
            pz_array[i]:=pz_array_n[i];

         x0:=x0n;
         dx_min:=dx_min_n;
         y0:=y0n;
         dy_min:=dy_min_n;
         z0:=z0n;
         dz_min:=dz_min_n;

         NodeCount_x:=New_NodeCount_x;
         NodeCount_y:=New_NodeCount_y;
         NodeCount_z:=New_NodeCount_z;


         setlength(SpaceChargeArray,New_NodeCount_x,New_NodeCount_y,New_NodeCount_z);

         Init_SpaceChargeArray(Voltage);

         for k := 0 to New_NodeCount_z-1 do
           for j := 0 to New_NodeCount_y-1 do
              for i := 0 to New_NodeCount_x-1 do
                   if SpaceChargeArray[i,j,k].Material<>mMetal then
                       begin
                          SpaceChargeArray[i,j,k].Phi:=tmp[i,j,k].Phi;
                          SpaceChargeArray[i,j,k].dPhi:=0;
                          SpaceChargeArray[i,j,k].n:=tmp[i,j,k].n;
                          SpaceChargeArray[i,j,k].p:=tmp[i,j,k].p;
                          SpaceChargeArray[i,j,k].dn:=0;
                          SpaceChargeArray[i,j,k].dp:=0;
                       end;

         setlength(tmp,0,0,0);
         output('refinement done.');
      end;

      procedure ExportIonizedDonorsToGnuPlotBinaryFile(Filename:String; y:integer);
      var ms:TMemoryStream;
          x,z:integer;
          px,py,pz:double;
          s,i:single;
          SC:integer;
          v:double;
      begin
        ms:=TMemoryStream.create;
        i:=NodeCount_x;
        ms.write(i,sizeof(i));

        for x:=0 to NodeCount_x-1 do
          begin
             GridToPhysical_Array(x,0,0,px,py,pz);
             s:=px;
             ms.write(s,sizeof(s));
          end;

        for z:=0 to NodeCount_z-1 do
          begin
              GridToPhysical_Array(0,0,z,px,py,pz);
              s:=pz;
              ms.write(s,sizeof(s));
              for x:=0 to NodeCount_x-1 do
                begin
                 if SpaceChargeArray[x,y,z].Material=mSC then
                    begin
                      SC:=SpaceChargeArray[x,y,z].SCIndex;
                      v:=ND_ionized(Semiconductors[SC].E_F, SpaceChargeArray[x,y,z].Phi+Semiconductors[SC].E_offset,SC);
                      s:=v;
                    end
                 else
                    s:=0;
                 ms.write(s,sizeof(s));
                end;
          end;

        ms.savetofile(filename);
        freeandnil(ms);
      end;

      // Saving/ Export procedure to store binary and ascii results
      procedure ExportIonizedAcceptorsToGnuPlotBinaryFile(Filename:String; y:integer);
      var ms:TMemoryStream;
          x,z:integer;
          px,py,pz:double;
          s,i:single;
          SC:integer;
          v:double;
      begin
        ms:=TMemoryStream.create;
        i:=NodeCount_x;
        ms.write(i,sizeof(i));

        for x:=0 to NodeCount_x-1 do
          begin
             GridToPhysical_Array(x,0,0,px,py,pz);
             s:=px;
             ms.write(s,sizeof(s));
          end;

        for z:=0 to NodeCount_z-1 do
          begin
              GridToPhysical_Array(0,0,z,px,py,pz);
              s:=pz;
              ms.write(s,sizeof(s));
              for x:=0 to NodeCount_x-1 do
                begin
                 if SpaceChargeArray[x,y,z].Material=mSC then
                    begin
                      SC:=SpaceChargeArray[x,y,z].SCIndex;
                      v:=NA_ionized(Semiconductors[SC].E_F, SpaceChargeArray[x,y,z].Phi+Semiconductors[SC].E_offset,SC);
                      s:=v;
                    end
                 else
                    s:=0;
                 ms.write(s,sizeof(s));
                end;
          end;

        ms.savetofile(filename);
        freeandnil(ms);
      end;

      // Saving/ Export procedure to store binary and ascii results
      procedure ExportPhiToGnuPlotBinaryFile(Filename:String; y:integer; SwapPhi:boolean);
      var ms:TMemoryStream;
          x,z:integer;
          px,py,pz:double;
          s,i:single;
          f:double;
      begin
        if SwapPhi then f:=-1 else f:=1;

        ms:=TMemoryStream.create;
        i:=NodeCount_x;
        ms.write(i,sizeof(i));

        for x:=0 to NodeCount_x-1 do
          begin
             GridToPhysical_Array(x,0,0,px,py,pz);
             s:=px;
             ms.write(s,sizeof(s));
          end;

        for z:=0 to NodeCount_z-1 do
          begin
              GridToPhysical_Array(0,0,z,px,py,pz);
              s:=pz;
              ms.write(s,sizeof(s));
              for x:=0 to NodeCount_x-1 do
                begin
                 s:=SpaceChargeArray[x,y,z].Phi*f;
                 ms.write(s,sizeof(s));
                end;
          end;

        ms.savetofile(filename);
        freeandnil(ms);
      end;

      // Saving/ Export procedure to store binary and ascii results
      procedure ExportPhiToAsciiFile(Filename:string;y:integer);
      var tf:Textfile;
          aline:string;
          x,z:integer;
          px,py,pz:double;
      begin
            assignfile(tf,filename);
            rewrite(tf);

            //Koordinatenachsen als erste Zeile
            aline:=inttostr(NodeCount_z)+',';
            for x:=0 to NodeCount_x-1 do
             begin
                GridToPhysical_Array(x,0,0,px,py,pz);
                aline:=aline+floattostr(px)+',';
             end;
            setlength(aline,length(aline)-1);
            writeln(tf,aline);

            for z := 0 to NodeCount_z-1 do
             begin
              GridToPhysical_Array(0,0,z,px,py,pz);
              aline:=floattostr(pz);
              for x := 0 to NodeCount_x-1 do
                begin
                  if length(aline)>0 then
                       aline:=aline+',';
                  aline:=aline+floattostr(spacechargearray[x,y,z].phi);
                end;
              writeln(tf,aline);
             end;
             closefile(tf);
      end;

      // Saving/ Export procedure to store binary and ascii results
      procedure ExportElectronDensityToGnuPlotBinaryFile(Filename:String; y:integer);
      var ms:TMemoryStream;
          x,z:integer;
          px,py,pz:double;
          s,i:single;
          d:double;
      begin

        ms:=TMemoryStream.create;
        i:=NodeCount_y;
        ms.write(i,sizeof(i));

        for x:=0 to NodeCount_x-1 do
          begin
             GridToPhysical_Array(x,0,0,px,py,pz);
             s:=px;
             ms.write(s,sizeof(s));
          end;

        for z:=0 to NodeCount_z-1 do
          begin
              GridToPhysical_Array(0,0,z,px,py,pz);
              s:=pz;
              ms.write(s,sizeof(s));
              for x:=0 to NodeCount_x-1 do
                begin
                 d:=SpaceChargeArray[x,y,z].n;
                 s:=d;
                 ms.write(s,sizeof(s));
                end;
          end;

        ms.savetofile(filename);
        freeandnil(ms);

      end;

      procedure ExportPhiGammaNToGnuPlotBinaryFile(Filename:String; y:integer; SwapPhi:boolean);
      var ms:TMemoryStream;
          x,z,f:integer;
          px,py,pz:double;
          s,i:single;
          d:double;
      begin
        if SwapPhi then f:=-1 else f:=1;

        ms:=TMemoryStream.create;
        i:=NodeCount_y;
        ms.write(i,sizeof(i));

        for x:=0 to NodeCount_x-1 do
          begin
             GridToPhysical_Array(x,0,0,px,py,pz);
             s:=px;
             ms.write(s,sizeof(s));
          end;

        for z:=(NodeCount_z div 2) to NodeCount_z-1 do
          begin
              GridToPhysical_Array(0,0,z,px,py,pz);
              s:=pz;
              ms.write(s,sizeof(s));
              for x:=0 to NodeCount_x-1 do
                begin
                 d:=SpaceChargeArray[x,y,z].Phi_Gamma_n*f;
                 s:=d;
                 ms.write(s,sizeof(s));
                end;
          end;

        ms.savetofile(filename);
        freeandnil(ms);

      end;

      procedure ExportPhiGammaPToGnuPlotBinaryFile(Filename:String; y:integer; SwapPhi:boolean);
      var ms:TMemoryStream;
          x,z,f:integer;
          px,py,pz:double;
          s,i:single;
          d:double;
      begin
        if SwapPhi then f:=-1 else f:=1;

        ms:=TMemoryStream.create;
        i:=NodeCount_y;
        ms.write(i,sizeof(i));

        for x:=0 to NodeCount_x-1 do
          begin
             GridToPhysical_Array(x,0,0,px,py,pz);
             s:=px;
             ms.write(s,sizeof(s));
          end;

        for z:=(NodeCount_z div 2) to NodeCount_z-1 do
          begin
              GridToPhysical_Array(0,0,z,px,py,pz);
              s:=pz;
              ms.write(s,sizeof(s));
              for x:=0 to NodeCount_x-1 do
                begin
                 d:=SpaceChargeArray[x,y,z].Phi_Gamma_p*f;
                 s:=d;
                 ms.write(s,sizeof(s));
                end;
          end;

        ms.savetofile(filename);
        freeandnil(ms);

      end;


      // Saving/ Export procedure to store binary and ascii results
      procedure ExportHoleDensityToGnuPlotBinaryFile(Filename:String; y:integer);
      var ms:TMemoryStream;
          x,z:integer;
          px,py,pz:double;
          s,i:single;
          d:double;
      begin

        ms:=TMemoryStream.create;
        i:=NodeCount_y;
        ms.write(i,sizeof(i));

        for x:=0 to NodeCount_x-1 do
          begin
             GridToPhysical_Array(x,0,0,px,py,pz);
             s:=px;
             ms.write(s,sizeof(s));
          end;

        for z:=0 to NodeCount_z-1 do
          begin
              GridToPhysical_Array(0,0,z,px,py,pz);
              s:=pz;
              ms.write(s,sizeof(s));
              for x:=0 to NodeCount_x-1 do
                begin
                 d:=SpaceChargeArray[x,y,z].p;
                 s:=d;
                 ms.write(s,sizeof(s));
                end;
          end;

        ms.savetofile(filename);
        freeandnil(ms);

      end;

      // Saving/ Export procedure to store binary and ascii results
      procedure SavePhiCentralToFile(Phi_central:array of tSpaceChargePoint; filename:string);
      var tf2:textfile;
          SC,i:integer;
      begin
             SC:=GetSemiconductorBelowTipApex;
             output('\nSaving central phi(z) to file: '+format(BB_z_filename,[Voltage])+'.');
             assignfile(tf2,format(BB_z_filename,[Voltage]));
             rewrite(tf2);
             writeln(tf2,'# z[nm], E_F[V], E_F+eV[eV], E_V[eV], E_C[eV]');
             for i := 0 to length(Phi_central)-1 do
               if Phi_central[i].Material=mSC then
                 with Semiconductors[SC] do
                  writeln(tf2,floattostr(pz_array[i])+', '+floattostr(E_F)+', '+floattostr(E_F+Voltage)+', '+floattostr(Phi_central[i].Phi)+', '+floattostr(Phi_central[i].Phi+E_g));
             closefile(tf2);
      end;

      // Saving/ Export procedure to store binary and ascii results
      procedure Save_n_and_p_ToFile(Phi_central:array of tSpaceChargePoint; filename:string);
      var tf2:textfile;
          i:integer;
      begin
             output('\nSaving central n(z) and p(z) to file: '+filename+'.');
             assignfile(tf2,filename);
             rewrite(tf2);
             writeln(tf2,'# z[nm], n[nm^-3], p[nm^-3]');
             for i := 0 to length(Phi_central)-1 do
               if Phi_central[i].Material=mSC then
                 writeln(tf2,floattostr(pz_array[i])+', '+floattostr(Phi_central[i].n)+', '+floattostr(Phi_central[i].p));
             closefile(tf2);
      end;

      // Saving/ Export procedure to store binary and ascii results
      procedure SaveIterimSolution(k:integer);
      var Dir:string;
      begin
        if (Interim_solution_Intervall>0) and (k mod Interim_solution_Intervall = 0)  then
                 begin
                  Dir:=IncludeTrailingPathDelimiter(format(Interim_solution_Path,[Voltage]));
                  if (k=0) and (not DirectoryExists(Dir)) then
                     mkdir(Dir);
                  ExportElectronDensityToGnuPlotBinaryFile(format(Dir+Interim_solution_n_filename+'.bin',[k]),Interim_solution_x);
                  ExportHoleDensityToGnuPlotBinaryFile(format(Dir+Interim_solution_p_filename+'.bin',[k]),Interim_solution_x);
                  ExportPhiToGnuPlotBinaryFile(format(Dir+Interim_solution_phi_filename+'.bin',[k]),Interim_solution_x,true);
                  ExportIonizedDonorsToGnuPlotBinaryFile(format(Dir+Interim_solution_ND_filename+'.bin',[k]),Interim_solution_x);
                  ExportIonizedAcceptorsToGnuPlotBinaryFile(format(Dir+Interim_solution_NA_filename+'.bin',[k]),Interim_solution_x);
                  ExportPhiGammaNToGnuPlotBinaryFile(format(Dir+Interim_solution_phi_quantum_n_filename+'.bin',[k]),Interim_solution_x,true);
                  ExportPhiGammaPToGnuPlotBinaryFile(format(Dir+Interim_solution_phi_quantum_p_filename+'.bin',[k]),Interim_solution_x,true);
                 end;
      end;



      // Swaping the electrostatic potential at each position x,y,z.
      // This is necessary before storing the results to a file, since Feenstra's
      // definition of the electrostatic potential (which is commonly used by
      // STM scientists) exhibits an opposing sign as compared to the commonly
      // known convention that is used in Selberherr's approach.
      // This can be directly seen when comparing the right hand side
      // of the Poisson quation that is used by Selberherr and Feenstra.
      // (In the case of Feenstra, there is a minus sign missing at the charge-
      // density.)
      procedure SwapPhi;
      var i,j,k:integer;
      begin
         for i := 0 to NodeCount_x-1 do
           for j := 0 to NodeCount_y-1 do
             for k := 0 to NodeCount_z-1 do
                 SpaceChargeArray[i,j,k].Phi:=-SpaceChargeArray[i,j,k].Phi;
      end;

       procedure ProgressUpdate(k,m:integer;CheckPoint_Phi,relChange,T:double); stdcall;
       begin
           output(format('k: %d |m: %d |phi_s: %.5fV |rel. change: %.3e |time: %.3fs',[k,m,Checkpoint_Phi,relChange,T]));
           SaveIterimSolution(k);
       end;


// Start of Pascal main program routine.
var
        i,j,x,y:integer;
        phi_central: array of TSpaceChargePoint;
        ivb,icb,iges:double;
        NumVoltage: integer;
        tf,tf2,tf3:textfile;
        Voltage_min,Voltage_max:double;
        x_central,y_central:integer;
        phi_s_central:double;
        listr:string;
        n_max,p_max,phi_max:double;
        SCbelowTip:integer;
        PhiScale_SC:double;

begin

        try
          output('Welcome to P_SpaceCharge!\n');
          output('=========================================================');
          output('Program version = ' + ProgramVersion);
          output('Author          = ' + 'Michael Schnedler');
          output('Target OS       = ' + TargetOS);
          output('=========================================================');
          sleep(500);
          output('\n');

          SetPrecisionMode(pmDouble);
          FormatSettings.Decimalseparator := '.';
          Ignore_1D_estimation := false;

          Quiet:= (UpperCase(ParamStr(1))='-Q');

          if not ReadParameterFile('parameters.txt') then
          begin
            output('\nPress return to exit program.');
            if not Quiet then readln;
            exit;
          end;




          Initialize_Semiconductors;
          ProgressNotification:=ProgressUpdate;


          NumVoltage:=length(V_List);

          Voltage_min:=1000;
          Voltage_max:=-1000;
          for j := 0 to NumVoltage-1 do
            begin
              if V_List[j]<Voltage_min then Voltage_min:=V_List[j];
              if V_List[j]>Voltage_max then Voltage_max:=V_List[j];
            end;

          if NumVoltage=0 then
          begin
            output('\nNo voltage point to calculate. Program Done.\nPress return to exit program.');
            if not Quiet then readln;
            exit;
          end;

          if I_V_output then
          begin
            assignfile(tf,I_V_filename);
            rewrite(tf);
            writeln(tf,'# V [V], I_vb [A], I_cb [A], I_ges [A]');
            closefile(tf);
          end;

          if BB_V_output then
          begin
             assignfile(tf2,BB_V_filename);
             rewrite(tf2);
             writeln(tf2,'#V[V], E_V[eV], E_C[eV], E_F + eV[eV]');
             closefile(tf2);
          end;

          if EFQ_V_output then
          begin
             assignfile(tf3,EFQ_V_Filename);
             rewrite(tf3);
             writeln(tf3,'#V[V], E_FQ_C [eV], E_FQ_V [eV]');
             closefile(tf3);
          end;


          for j := 0 to NumVoltage-1 do
            begin
               Voltage:=V_List[j];
               output(format('\nCalculation of the 3D potential solution (BIAS = %.5e V):\n',[voltage]));
               Init(Voltage);

               SCbelowTip:=GetSemiconductorBelowTipApex;

               for i := 0 to refinement_steps do
                begin
                  Iterate;
                  if i<refinement_steps then
                     RefineGrid(NodeCount_x*2-1,NodeCount_y*2-1,NodeCount_z*2-1);
                end;

               x_central:=NodeCount_x div 2;
               y_central:=NodeCount_y div 2;


               if Semiconductors[SCbelowTip].light_on then
                 begin
                    E_FQ_C:= GetE_Fqn2(SCbelowTip,Phi(x_central,y_central,NodeCount_z div 2)+Semiconductors[SCbelowTip].E_offset,Get_n(x_central,y_central,NodeCount_z div 2));
                    E_FQ_V:= GetE_Fqp2(SCbelowTip,Phi(x_central,y_central,NodeCount_z div 2)+Semiconductors[SCbelowTip].E_offset,Get_p(x_central,y_central,NodeCount_z div 2));
                 end
               else
                 begin
                    E_FQ_C:=Semiconductors[SCbelowTip].E_F;
                    E_FQ_V:=Semiconductors[SCbelowTip].E_F;
                 end;

               if EFQ_V_output then
               begin
                  output('\nWriting to file: '+EFQ_V_Filename);
                  assignfile(tf3,EFQ_V_Filename);
                  append(tf3);
                  writeln(tf3,floattostr(Voltage)+','+ floattostr(E_FQ_C)+','+floattostr(E_FQ_V));
                  closefile(tf3);
               end;

               SwapPhi;

               setlength(Phi_central, NodeCount_z);

               for i := 0 to NodeCount_z-1 do
                 begin
                   Phi_central[i]:=SpaceChargeArray[x_central,y_central, i];
                   Phi_central[i].Phi:=Phi_central[i].Phi-Semiconductors[SCbelowTip].E_offset;
                                                         // minus E_offset due to SwapPhi
                   if (i>0) and (SpaceChargeArray[x_central,y_central,i].Material=mSC)
                       and (SpaceChargeArray[x_central,y_central,i-1].Material=mVac) then
                            phi_s_central:=Phi_central[i].Phi;
                 end;

               iges:=I_Tunnel(Voltage,Phi_Central,pz_array,ivb,icb,Semiconductors[GetSemiconductorBelowTipApex]);

              if I_V_output then
                begin
                  output('\nWriting to file: '+I_V_filename);
                  assignfile(tf,I_V_filename);
                  append(tf);
                  writeln(tf, floattostr(voltage)
                              +','+floattostr(ivb)+','+floattostr(icb)
                              +','+floattostr(iges));
                  closefile(tf);
                end;

               if BB_z_output then
                 SavePhiCentralToFile(Phi_central,format(BB_z_filename,[Voltage]));

               if Rho_central_output then
                  Save_n_and_p_ToFile(Phi_central,format(Rho_central_filename,[Voltage]));

               if BB_V_output then
               begin
                assignfile(tf2,BB_V_filename);
                append(tf2);
                output('\nWriting to file: '+BB_V_Filename);
                writeln(tf2,floattostr(Voltage)+', '+
                            floattostr(Phi_s_central)+', '+
                            floattostr(Phi_s_central+Semiconductors[SCbelowTip].E_g)+', '+
                            floattostr(Semiconductors[SCbelowTip].E_F+Voltage));
                closefile(tf2);
               end;

               if Phi3D_output then
                begin
                 output('\nWriting Phi3D Arrays.');
                 for y := 0 to NodeCount_y-1 do
                  begin
                    ExportPhiToGnuPlotBinaryFile(format(Phi3D_filename+'.bin',[Voltage,y]),y,false);
                    ExportPhiToAsciiFile(format(Phi3D_filename+'.txt',[Voltage,y]),y);
                  end;
                end;

               if Rho3D_output then
                begin
                 output('\nWriting Rho3D Arrays.');
                 for y := 0 to NodeCount_y-1 do
                  begin
                    ExportElectronDensityToGnuPlotBinaryFile(format(Rho3D_n_filename+'.bin',[Voltage,y]),y);
                    ExportHoleDensityToGnuPlotBinaryFile(format(Rho3D_p_filename+'.bin',[Voltage,y]),y);
                  end;
                end;


            end;

            output('\n\nProgram done.');
            if not Quiet then readln;
        except
          on E: Exception do
            begin
              Writeln(E.ClassName, ': ', E.Message);
              if not Quiet then readln;
            end;
        end;
end.
