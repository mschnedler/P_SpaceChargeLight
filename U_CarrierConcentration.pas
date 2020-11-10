unit U_CarrierConcentration;

interface
uses U_Constants, U_Common, SysUtils, Math;

  function  F_j(j:double;eta:double):double;
  procedure Init_FjList(Phi_Start,Phi_End:double; Count: integer);
  function  GetFromFjList(Phi:double):double;
  function  n(E_F,Phi:double;inv:boolean; SC:integer):double;
  function  p(E_F,Phi:double; inv:boolean; SC:integer):double;
  function  ND_ionized(E_F,Phi:double; SC:integer):double;
  function  dND_ionized_dPhi(E_F,Phi:double; SC:integer):double;
  function  NA_ionized(E_F,Phi:double; SC:integer):double;
  function  dNA_ionized_dPhi(E_F,Phi:double; SC:integer):double;
  function  rho(E_F:double; Phi:double; invc, invv:boolean; SC:integer):double;
  function  rho_GS(x: double; Parameters:TFunctionParameters): double;
  function  Find_EF(SC:integer):double;

  procedure Find_EF2(E_F_min,E_F_max,Precision:double; SC:integer);

  procedure Init_SurfRhoList(Count:integer; SC:integer);
  function  GetRhoSurf(Phi:double;var rho:double; SC:integer):boolean;
  function  drho_surf_dPhi(Phi:double; SC:integer):double;

  function GetE_Fqn2(SC:integer; Phi:double; n:double):double;
  function GetE_Fqp2(SC:integer; Phi:double; p:double):double;

implementation


       // Bibliography
       //  [1] Seiwatz and Green, J. Appl. Phys. 29 (7), 1958


       // *************************************************************************
       // ***  Carrier concentrations in parabolic band approx.                 ***
       // *************************************************************************

       // Helping function for Fj(eta), see below and [1]:
       function F_j_int(x:double; fp:TFunctionParameters):double;
       begin
         result:= power(x,fp[0])/(1+exp(x-fp[1]));
       end;

       // Definition of the Fermi-Dirac Integral Fj(eta), see [1]:
       function F_j(j:double;eta:double):double;
       var fp:TFunctionParameters;
       begin
           setlength(fp,2);
           fp[0]:=j;
           fp[1]:=eta;
           if eta > 40 then
             result:=power(eta,(j+1))/(j+1)
           else if eta < -8 then
             begin
               result:=sqrt(Pi)*exp(eta)/2;
               if abs(j-3/2)<0.0001 then result:=result*3/2;
             end
           else
             result:=integrate(F_j_int,fp,0,20+eta,Integration_Accuracy);  //IMPORTANT: IF THERE ARE
         end;                                               //CONVERGENCE PROBLEMS IN THE
                                                            //POTENTIALDETERMINATION: SET
                                                            //THE UPPER LIMIT 20+eta TO HIGHER
                                                            //VALUES BECAUSE THIS INTEGRAL
                                                            //USUALLY GOES TO INFINITY!

      procedure Init_FjList(Phi_Start,Phi_End:double; Count: integer);
      var i,j:integer;
      begin
        if t=0 then exit;
        output(format('\nPreparing FjList (from Phi= %.5e V to Phi= %.5e V)',[Phi_Start,Phi_End]));
        setlength(FjList.Value,count);
        FjList.Phi_Start:=Phi_Start;
        FjList.Phi_End:=Phi_End;
        FjList.Count:=Count;
        for j := 0 to Count-1 do
             FjList.Value[j]:= F_j(0.5,(Phi_Start+(Phi_End-Phi_Start)*j/(Count-1))/(k*T));
        output('Preparing FjList...done');
      end;

      function GetFromFjList(Phi:double):double;
      var n1,n2:integer;
          n:double;
          test:double;
      begin
        n:=(FjList.Count-1)*(Phi-FjList.Phi_Start)/(FjList.Phi_End-FjList.Phi_Start);
        n1:=floor(n);
        n2:=ceil(n);
        if (n1>=0) and (n2<FjList.Count) then
           begin
             result:=FjList.Value[n1]*(n2-n)+FjList.Value[n2]*(n-n1);
           end
        else
          result:=F_j(0.5,Phi/(k*T));
      end;

       // Definition of n(E), the electrons in the conduction band
       // n returns the electron concentration in [1/nm^3]
       // Phi is given in Selberherr's sign convention
       function n(E_F,Phi:double;inv:boolean; SC:integer):double;
       var c:double;
       begin
           //*************************************************************************************************
           //* c:=2/sqrt(Pi)*2*power(me/(2*Pi*h_bar*h_bar),1.5)*1e-27*power(eVToJ,1.5); // in eV^-1.5*nm^-3  *
           //* In order to optimize computation time, we store the result of                                 *
           //*************************************************************************************************
           c:=6.81232110559823;
           if not inv then Phi:=0;
           if T=0 then
           begin
             if (E_F-Semiconductors[SC].E_g+Phi<0) then
               result:=0
             else
               result:= 2*C/3*power(Semiconductors[SC].m_cb*(E_F-Semiconductors[SC].E_g+Phi),1.5);
             exit;
           end;

           result:=C*power(Semiconductors[SC].m_cb*kb*T/eVToJ,1.5)*GetFromFjList(E_F-Semiconductors[SC].E_g+Phi);
          // c:=2*me*kb*T/(h*h)*1e-18; //[1/nm^2]
           //if not inv then result:=0
           // else

           //result:=4*Pi*power(Semiconductors[SC].m_cb*c,3/2)*GetFromFjList(E_F-Semiconductors[SC].E_g+Phi);
           //result:=4*Pi*power(Semiconductors[SC].m_cb*c,3/2)*F_j(0.5,(E_F-Semiconductors[SC].E_g+Phi)/(k*T));
       end;


       // Definition of p(E), the holes in the valence band, see [1]:
       // p returns the hole concentration in [1/nm^3]
       // Phi is given in Selberherr's sign convention
       function p(E_F,Phi:double; inv:boolean; SC:integer):double;
       var c:double;
       begin
           //*************************************************************************************************
           //* c:=2/sqrt(Pi)*2*power(me/(2*Pi*h_bar*h_bar),1.5)*1e-27*power(eVToJ,1.5); // in eV^-1.5*nm^-3  *
           //* In order to optimize computation time, we store the result of                                 *
           //*************************************************************************************************
           c:=6.81232110559823;
           if not inv then Phi:=0;
           if T=0 then
           begin
             if (-Phi-E_F<=0) then
               result:=0
             else
               result:= 2*C/3*power(Semiconductors[SC].m_vb*(-Phi-E_F),1.5);
             exit;
           end;
           result:=C*power(Semiconductors[SC].m_vb*kb*T/eVToJ,1.5)*GetFromFjList(-Phi-E_F);
        //   c:=2*me*kb*T/(h*h)*1e-18; //[1/nm^2]
        //   result:=4*Pi*power(Semiconductors[SC].m_vb*c,3/2)*GetFromFjList(-Phi-E_F);
//           result:=4*Pi*power(Semiconductors[SC].m_vb*c,3/2)*F_j(0.5,(-Phi-E_F)/(k*T));
       end;

       // Definition of ionized donors, see[1]:
       // Returns the concentration of ionized donors in [1/nm^3]
       // Phi is given in Selberherr's sign convention
       function ND_ionized(E_F,Phi:double; SC:integer):double;
       var Expo:double;
       begin
         Expo:=(E_F - (Semiconductors[SC].E_D-Phi));
         if T=0 then
         begin
            if Expo<=0 then
              result:= Semiconductors[SC].C_ND
            else
              result:=0;
            exit;
         end;
         Expo:=Expo/(k*T);
         if Expo<-40 then
            result:=Semiconductors[SC].C_ND
         else if Expo>40 then
            result:=0
         else
            result:= Semiconductors[SC].C_ND * 1/(1+2*exp(Expo));

       end;

       // Definition of the derivative of the ionized donors with respect to Phi
       // dND_ionized/dPhi in [1/(nm^3*eV)]
       // Phi is given in Selberherr's sign convention
       function dND_ionized_dPhi(E_F,Phi:double; SC:integer):double;
       var Expo:double;
       begin
         Expo:=E_F-Semiconductors[SC].E_D+Phi;
         if T=0 then
           result:=0
         else
           begin
             Expo:=Expo/(k*T);
             if (Expo<-40) or (Expo>40) then
               result:=0
             else
               result:= - Semiconductors[SC].C_ND * 2 *exp(Expo)/(k*T*power(2*exp(Expo)+1,2));
           end;
       end;

       // Definition of ionized acceptors, see[1]:
       // Returns the concentration of ionized acceptors in [1/nm^3]
       // Phi is given in Selberherr's sign convention
       // The factor of 4 in front of the exponential function
       // takes into account two-fold degenerate valence band, which is a valid
       // assumption for most commonly used semiconductors.
       function NA_ionized(E_F,Phi:double; SC:integer):double;
       var Expo:double;
       begin
         Expo:=((Semiconductors[SC].E_A-Phi) - E_F);
         if T=0 then
         begin
          if Expo<=0 then
            result:=Semiconductors[SC].C_NA
          else
            result:=0;
          exit;
         end;
         Expo:=Expo/(k*T);
         if Expo<-40 then
           result:=Semiconductors[SC].C_NA
         else if Expo>40 then
           result:=0
         else
           result:= Semiconductors[SC].C_NA * 1/(1+4*exp(Expo));

       end;

       // Definition of the derivative of the ionized acceptors with respect to Phi
       // dNA_ionized/dPhi in [1/(nm^3*eV)]
       // Phi is given in Selberherr's sign convention
       function dNA_ionized_dPhi(E_F,Phi:double; SC:integer):double;
       var Expo:double;
       begin
         Expo:=Semiconductors[SC].E_A-E_F-Phi;
         if T=0 then
           result:=0
         else
           begin
             Expo:=Expo/(k*T);
             if (Expo<-40) or (Expo>40) then
               result:=0
             else
               result:=Semiconductors[SC].C_NA * 2 *exp(Expo)/(k*T*power(2*exp(Expo)+1,2));
           end;
       end;

       // Definition of the whole charge carrier concentration, see[1]
       // Returns the charge carrier concentration in [1/nm^3]
       function rho(E_F:double; Phi:double; invc, invv:boolean; SC:integer):double;
       var n0,p0:double;
       begin
         n0:=n(E_F,Phi,invc,SC);
         p0:=p(E_F,Phi,invv,SC);

         result:=   (ND_ionized(E_F,Phi,SC)-
                     NA_ionized(E_F,Phi,SC)+
                     p0-
                     n0);
       end;


       //Wrapper of function rho for the determination of E_F (Golden section
       //iteration). Do not use elsewhere!
       function  rho_GS(x: double; Parameters:TFunctionParameters): double;
       var SC:integer;
       begin
          SC:=trunc(Parameters[0]);
          result:=abs(rho(x, 0, true, true, SC));
       end;

       // Determination of the Fermi-Energy E_F [eV]
       // Numerically solves the charge neutrality condition in thermal equilibrium
       // rho(E_F) = e*(N_D - N_A + p - n) = 0
       // Algorithm used for determination of the zero-point is Bisection Iteration
       // E_F_min and E_F_max define the interval for the Bisection algorithm
       function Find_EF(SC:integer):double;
       var i:integer;
           E_F_min,E_F_max,dE:double;
           Param:TFunctionParameters;
           r_max,r_min, curr_r: double;
           Num_iter:integer;
           tf:textfile;
       begin
         //Intrinsic case:
         if (Semiconductors[SC].N_D=0) and (Semiconductors[SC].N_A=0) then
         begin
           result:= Semiconductors[SC].E_G/2 +
                    0.75*kb*T*ln(Semiconductors[SC].m_vb/Semiconductors[SC].m_cb)/eVToJ;
           output(format('Fermi level found at E_F - E_V = %.5e eV.\n',[result]));
           output(format('Electron concentration = %.5e nm^-3',[n(result,0,true,SC)]));
           output(format('Hole concentration = %.5e nm^-3',[p(result,0,true,SC)]));
           output(format('Ionized donor concentration = %.5e nm^-3',[ND_ionized(result,0,SC)]));
           output(format('Ionized acceptor concentration = %.5e nm^-3\n',[NA_ionized(result,0,SC)]));
           exit;
         end;

         setlength(Param,1);
         Param[0]:=SC;
         //T=0K:
         if T=0 then
           begin
             if Semiconductors[SC].N_A>Semiconductors[SC].N_D then
               begin
                 E_F_min:=-1;
                 E_F_max:=0;
               end
             else if Semiconductors[SC].N_A<Semiconductors[SC].N_D then
               begin
                 E_F_min:=Semiconductors[SC].E_G;
                 E_F_max:=Semiconductors[SC].E_G+1;
               end
             else
               begin //Semiconductors[SC].N_A=Semiconductors[SC].N_D
                 result:=(Semiconductors[SC].E_D+Semiconductors[SC].E_A)/2;
                 output(format('Fermi level found at E_F - E_V = %.5e eV.\n',[result]));
                 output(format('Electron concentration = %.5e nm^-3',[n(result,0,true,SC)]));
                 output(format('Hole concentration = %.5e nm^-3',[p(result,0,true,SC)]));
                 output(format('Ionized donor concentration = %.5e nm^-3',[ND_ionized(result,0,SC)]));
                 output(format('Ionized acceptor concentration = %.5e nm^-3\n',[NA_ionized(result,0,SC)]));
                 exit;
               end;
           end
         else
           begin //T>0K
             //Estimation of E_F:
             E_F_max:=-0.3;
             r_max:=abs(rho(E_F_max,0,true,true,SC));
             dE:=(Semiconductors[SC].E_G+0.6)/5000;
             for i := 0 to 4999 do
               begin
                 E_F_max:=E_F_max+dE;
                 curr_r:=abs(rho(E_F_max,0,true,true,SC));
                 if curr_r>r_max then break;
                 r_max:=curr_r;
               end;
             E_F_min:=E_F_max-2*dE;
           end;

        { assignfile(tf,'test.txt');
         rewrite(tf);
         for i := 0 to 999 do
          begin
            curr_r:=rho(E_F_min+(E_F_max-E_F_min)/999*i,0,true,true,SC);
            writeln(tf,floattostr(E_F_min+(E_F_max-E_F_min)/999*i)+', '+floattostr(curr_r)+', '+floattostr(NA_ionized(E_F_min+(E_F_max-E_F_min)/999*i,0,SC))+', '+floattostr(p(E_F_min+(E_F_max-E_F_min)/999*i,0,true, SC)));
          end;
         closefile(tf);}
         output('\nE_F for semiconductor '+inttostr(SC)+':');
         output(Format('Search interval: [%.5e eV, %.5e eV]',[E_F_min,E_F_max]));
         Num_Iter:= GoldenSectionInteration(rho_GS,param, E_F_min, E_F_max,Precision);
         result:=(E_F_min+E_F_max)/2;
         output(format('Fermi level found at E_F - E_V = %.5e eV with %d iterations.\n',[result,Num_Iter]));
         output(format('Electron concentration = %.5e nm^-3',[n(result,0,true,SC)]));
         output(format('Hole concentration = %.5e nm^-3',[p(result,0,true,SC)]));
         output(format('Ionized donor concentration = %.5e nm^-3',[ND_ionized(result,0,SC)]));
         output(format('Ionized acceptor concentration = %.5e nm^-3\n',[NA_ionized(result,0,SC)]));

       end;


       procedure Find_EF2(E_F_min,E_F_max,Precision:double; SC:integer);
       var i:integer;
           r_min,r_max,curr_r:double;
       begin
         i:=0;
         output('\nBisection iteration of E_F for semiconductor '+inttostr(SC)+':\n');
         r_min:=rho(E_F_min,0,true,true,SC);
         r_max:=rho(E_F_max,0,true,true,SC);
         if sign(r_min)<>sign(r_max) then
          while (abs(E_F_max-E_F_min)>Precision) and (i<newton_steps_max) do
            begin
              curr_r:=rho(E_F_min+(E_F_max-E_F_min)/2,0,true,true,SC);
              if curr_r=0 then
                begin
                  // already found the correct result
                  // important for intrinsic sc
                  Semiconductors[SC].E_F:=(E_F_max+E_F_min)/2;
                  E_F_max:=Semiconductors[SC].E_F;
                  E_F_min:=Semiconductors[SC].E_F;
                  break;
                end;
              if sign(curr_r)=sign(r_min) then
                  E_F_min:=E_F_min+(E_F_max-E_F_min)/2
              else
                  E_F_max:=E_F_min+(E_F_max-E_F_min)/2;
              inc(i);
              r_min:=rho(E_F_min,0,true,true,SC);
              r_max:=rho(E_F_max,0,true,true,SC);
              Semiconductors[SC].E_F:=(E_F_max+E_F_min)/2;
              output(format('(%d) E_F - E_V = %.5e eV',[i,Semiconductors[SC].E_F]));
            end;
         if (i=newton_steps_max) and (abs(E_F_max-E_F_min)>Precision) then
            begin
               raise EMathError.Create(format('Bisection iteration does not converge:\nE_F not found within %d steps and precision = %.5e.\nTry to increase newton_steps_max or decrease precision!',[newton_steps_max,precision]));
            end
         else
            output(format('\nE_F found at (E_F - E_V) = %.5e eV \n(rho_min = %.5e C/m^3, rho_max = %.5e C/m^3)',[Semiconductors[SC].E_F,r_min,r_max]));
       end;



       // Definition of the Gaussian distribution (normal distribution)
       // for the calculation of surface states
       function GaussDist(E:double;p:TFunctionParameters):double;
       var E_SS:double;
       begin
         E_SS:=p[2];
         result:=exp(-0.5*power((E-E_SS)/p[1],2));
       end;

       // Integration over Gaussian distribution
       // Returns the surface charge density in [C/m^2]
       // Integration from charge neutrality level E_CNL to E_F
       // The position of E_SS (peak position of normal distribution) changes
       // with Phi. To get an array containing the surface charge density for
       // an arbitrary value, Phi is taken out of the integrand and into the
       // limits of the integration of the Gaussian distribution:
       // integral from E_CNL+Phi to E_F -> G(E-(E_SS+Phi)) dE =
       // integral from E_CNL to E_F-Phi -> G(E-E_SS) dE
       // Is used for initializing the array RhoSurfList, that stores the
       // surface charge density as a function of the potential Phi.
       // Without this precalculation, the poisson solver would
       // slow down significantly.
       // This verion of rho_surf is optimized in order fill an array. It uses
       // the previous entry of the array to calculate the next value.
       function rho_surf(Rho_prev,Phi_prev,Phi:double; SC:integer):double;
       var fp:TFunctionParameters;
           w:double;
           c:double;
       begin
           if Semiconductors[SC].surface_charge_density=0 then
            begin
              result:=0;
              exit;
            end;
           w:=Semiconductors[SC].FWHM/(2*sqrt(2*ln(2)));
           setlength(fp,3);
           fp[0]:=Phi;
           fp[1]:=w;
           fp[2]:=Semiconductors[SC].E_SS;
           c:=-e*Semiconductors[SC].surface_charge_density/sqrt(2*Pi*w*w);
           result:=rho_prev+c*integrate(GaussDist,fp,Semiconductors[SC].E_F-Phi_prev,Semiconductors[SC].E_F-Phi,1e-10);
       end;

       // Definition of the derivative of the surface charge density with respect to Phi
       // drho_surf/dPhi in [C/(m^2*eV)]
       function drho_surf_dPhi(Phi:double; SC:integer):double;
       var w:double;
       begin
        result:=0;
        if Semiconductors[SC].surface_charge_density=0 then exit;
        w:=Semiconductors[SC].FWHM/(2*sqrt(2*ln(2)));
        result:=exp(-0.5*power((Semiconductors[SC].E_F-Semiconductors[SC].E_SS-Phi)/w,2))*Semiconductors[SC].surface_charge_density/sqrt(2*Pi*w*w)*e;
       end;



      // Initializing the surface charge density array.
      // Uses the procedure rho_surf, as described above.
      // An interval of 10 sigma around the surface state distribution is used.
      // The initial value of Phi_prev is choosen such, that the lower
      // integration boundary of the integration in rho_surf is
      // set to E_CNL (see comments before rho_surf for further information).
      procedure Init_SurfRhoList(Count:integer; SC:integer);
      var Phi_Start,Phi_End:double;
          w:double;
          i: Integer;
          Phi_new,Phi_prev,Rho_prev:double;
      begin
         output('\nPreparing surface charge density of semiconductor '+inttostr(SC)+'.\n');
         w:=Semiconductors[SC].FWHM/(2*sqrt(2*ln(2)));
         Phi_Start:=Semiconductors[SC].E_F-Semiconductors[SC].E_SS-10*w;
         Phi_End:=Semiconductors[SC].E_F-Semiconductors[SC].E_SS+10*w;
         Semiconductors[SC].SurfRhoList.Count:=Count;
         Semiconductors[SC].SurfRhoList.Phi_Start:=Phi_Start;
         Semiconductors[SC].SurfRhoList.Phi_End:=Phi_End;
         setlength(Semiconductors[SC].SurfRhoList.Value,Count);

         Phi_prev:=Semiconductors[SC].E_F-Semiconductors[SC].E_CNL;
         Rho_prev:=0;

         for i := 0 to Count-1 do
           begin
             if Semiconductors[SC].surface_charge_density=0 then
                Semiconductors[SC].SurfRhoList.Value[i]:=0
              else
                begin
                  Phi_new:=Phi_Start+(Phi_End-Phi_Start)/(Count-1)*i;
                  Rho_prev:=rho_surf(Rho_prev,Phi_prev,Phi_new,SC);
                  Semiconductors[SC].SurfRhoList.Value[i]:=Rho_prev;
                  Phi_prev:=Phi_new;
                end;
           end;
      end;

      // Procedure for receiving the surface charge distribution for a given
      // Phi from the previously filled array. This procedure is used in the
      // latter finite-difference iteration to obtain the surface charge
      // densities fast enough.
      // When Phi exceeds the borders of the array, the last value of the array
      // is used as return value. If Phi is in between two array entries,
      // a linear approximation between the two nearest array entries is used.
      // The result is given in units of [C/m^2]
      function GetRhoSurf(Phi:double;var rho:double; SC:integer):boolean;
      var n1,n2:integer;
          n:double;
          E:double;
      begin
        result:=true;

        if (Semiconductors[SC].surface_charge_density=0) or
           (Semiconductors[SC].SurfRhoList.Phi_Start=Semiconductors[SC].SurfRhoList.Phi_End) then
         begin
           rho:=0;
           exit;
         end;

        E:=Phi;

        n:=(Semiconductors[SC].SurfRhoList.Count-1)*(E-Semiconductors[SC].SurfRhoList.Phi_Start)/
           (Semiconductors[SC].SurfRhoList.Phi_End-Semiconductors[SC].SurfRhoList.Phi_Start);
        n1:=floor(n);
        n2:=ceil(n);

        if (n1<0) then
             rho:=Semiconductors[SC].SurfRhoList.Value[0]
        else if (n2>=Semiconductors[SC].SurfRhoList.Count) then
             rho:=Semiconductors[SC].SurfRhoList.Value[Semiconductors[SC].SurfRhoList.Count-1]
        else
             rho:=Semiconductors[SC].SurfRhoList.Value[n1]+
                  (Semiconductors[SC].SurfRhoList.Value[n2]-Semiconductors[SC].SurfRhoList.Value[n1])*(n-int(n));

      end;


      //Returns the Quasi-Fermi level for electrons.
      //Phi is defined in Selberherr notation.
      //The quasi-Fermi level commonly shifts/moves with the potential phi.
      //Hence, the result depends on phi. In the result of this function, phi is given
      //in Feenstra's sign convention. This may be a little bit confusing. However,
      //all further calculations that depend on the quasi-Fermi levels are derived
      //in Feenstra's sign convention.
      function GetE_Fqn2(SC:integer; Phi:double; n:double):double;
       var x:double;
           x2,x3,x4:double;
      begin
           //Joyce-Dixon Approximation (see [6])
           //The Joyce-Dixon approximation is valid for all negative
           //values of x as it asymptotes to the low temperature/low density
           //limit but is only valid for positive values of x below about 5.
           x:=n/Semiconductors[SC].N_C+1e-20;

           if IsZero(x) then
             begin
               result:= Semiconductors[SC].E_F;//Semiconductors[SC].E_g-Phi-40*k*T;
               exit;
             end;

           if x<8 then
            begin
              x2:=x*x;
              x3:=x2*x;
              x4:=x2*x2;
              result:=(Semiconductors[SC].E_g-Phi)+k*T*(ln(x)+FDI_Coefficients[0]*x+FDI_Coefficients[1]*x2+FDI_Coefficients[2]*x3+FDI_Coefficients[3]*x4);

              //Smooth/linear transition between both approximations (otherwise, m-iteration may not reach required convergence critera):
              if x>7 then
                result:=result*(1-(x-7))+ ((Semiconductors[SC].E_g-Phi)+k*T*sqrt(power(3*sqrt(Pi)/4*x,4/3)-Pi*Pi/6)) * (x-7);
            end
           else
            result:=(Semiconductors[SC].E_g-Phi)+k*T*sqrt(power(3*sqrt(Pi)/4*x,4/3)-Pi*Pi/6);
      end;

      //Returns the quasi-Fermi level for holes.
      //Phi is defined in Selberherr notation.
      //The quasi-Fermi level commonly shifts/moves with the potential phi.
      //Hence, the result depends on phi. In the result of this function, phi is given
      //in Feenstra's sign convention. This may be a little bit confusing. However,
      //all further calculations that depend on the quasi-Fermi levels are derived
      //in Feenstra's sign convention.
      function GetE_Fqp2(SC:integer; Phi:double; p:double):double;
       var x:double;
           x2,x3,x4:double;
      begin
           // Joyce-Dixon Approximation (see [6])
           // The Joyce-Dixon approximation is valid for all negative
           // values of x as it asymptotes to the low temperature/low density
           // limit but is only valid for positive values of x below about 5.

           x:=p/Semiconductors[SC].N_V+1e-20;

           if IsZero(x) then
             begin
               result:= Semiconductors[SC].E_F;// -Phi+40*k*T;
               exit;
             end;


           if x<8 then
            begin
              x2:=x*x;
              x3:=x2*x;
              x4:=x2*x2;
              result:=(-Phi)-k*T*(ln(x)+FDI_Coefficients[0]*x+FDI_Coefficients[1]*x2+FDI_Coefficients[2]*x3+FDI_Coefficients[3]*x4);
              //Smooth/linear transition between both approximations (otherwise, m-iteration may not reach required convergence critera):
              if x>7 then
                result:=result*(1-(x-7))+ (-Phi-k*T*sqrt(power(3*sqrt(Pi)/4*x,4/3)-Pi*Pi/6)) * (x-7);
            end
           else
            result:=(-Phi)-k*T*sqrt(power(3*sqrt(Pi)/4*x,4/3)-Pi*Pi/6);
      end;


end.
