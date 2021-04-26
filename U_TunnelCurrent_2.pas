unit U_TunnelCurrent_2;
// Literature:
// [1] J. Vac. Sci. Technol. B 5 (923), 1897

// Calculates the tunnel current in the approximation of Feenstra[1].
// Latest modification: 31.07.2019
// Michael Schnedler

// History:
// 07.12.2016: Changed D_vac to D_vac_times_D_vb in all functions Theta_times_D_vac_dW* for the valence band
//             Changed D_vac to D_vac_times_D_cb in all functions Theta_times_D_vac_dW* for the conduction band
// 28.08.2019: Removed I_cb2 to I_cb5, I_vb2 to I_vb5 and implemented the integration over all VB and CB states
//             in I_VB and I_CB, respectively

{$IFDEF FPC}
  {$MODE DELPHI}
{$ENDIF}

interface

       uses U_Common,U_Constants, Math, SysUtils;

       function Sqrt_E_cb_minus_W(z:double; Parameters: TFunctionParameters):double;
       function D_CB(W:double; Parameters_In:TFunctionParameters):double;
       function Sqrt_W_minus_E_vb(z:double; Parameters: TFunctionParameters):double;
       function D_VB(W:double; Parameters_In:TFunctionParameters):double;
       function D_vac(W:double; Parameters_In:TFunctionParameters):double;
       function D_vac_times_D_cb(W:double;Parameters_In:TFunctionParameters):double;
       function D_vac_times_D_vb(W:double;Parameters_In:TFunctionParameters):double;

       function Theta_times_int_D_CB_D_Vac_dW(E:double;Parameters:TFunctionParameters):double;
       function Theta_times_int_D_VB_D_vac_dW(E:double;Parameters:TFunctionParameters):double;
       //function Theta_times_int_D_vac_dW(E:double;Parameters:TFunctionParameters):double;

       function I_CB(Voltage:double;Parameters:TFunctionParameters):double;

       function I_VB(Voltage:double;Parameters:TFunctionParameters):double;

       function I_tunnel(Voltage:double;  Phi_central:array of TSpaceChargePoint; pz_array:array of double; var ivb,icb:double; Semiconductor:TSemiconductor):double;


implementation
      var Area:double;
          SC:TSemiconductor;
          E_VB_Max,E_CB_Min: double;

//*********************************************
//*Calculation of the barriers for tunneling: *
//*********************************************

// 1. Tunneling into conduction band thru space charge region (upward band bending)

       // 1a. Helping function for the integration of the barrier-potential in space charge region of the conduction band
       // z: distance from SC surface
       // Parameters[0]: W
       // Parameters[1]: dz
       // Parameters[2]: num //Number of following array entries
       // Parameters[3..3+num]: phi-array for the integration as prepared in D_CB, Values in eV
       function Sqrt_E_cb_minus_W(z:double; Parameters: TFunctionParameters):double;
       var Re:double;
           num,i:integer;
           iz:double;
           dz:double;
           W:double;
           phi_x:double;
       begin
         W:=Parameters[0];
         dz:=Parameters[1];
         num:=trunc(Parameters[2]);

         //Linear interpolation of the potential for smooth integration
         iz:=z/dz;
         i:=trunc(iz);
         if i<num-1 then
           begin
             iz:=iz-i;
             phi_x:= Parameters[i+3]*(1-iz)+Parameters[(i+1)+3]*iz;
           end
         else
           phi_x:=Parameters[num-1+3];

         // if holes fill the valence band above E_F then the barrier for
         // electrons tunneling thru SCR (and valence band) into conduction band
         // is zero in valence band region.
         // if no holes are in valence band, barrier is given by conduction band edge

         Re := (SC.E_g+phi_x-W);
         if (SC.Inv_V) and (W<phi_x) and (W>SC.E_F) then
           result:=0
         else if Re>0 then
                result:=sqrt( Re )
              else
                result:=0;
       end;

       // 1b. calculates the transmission coefficient for tunneling thru the space charge region of the
       // conduction band. see Eq. 4, 5 and Section E in [2]
       // W: Energy in eV
       // Parameters_In[0]:  Voltage in V
       // Parameters_In[1]:  z_max_sc
       // Parameters_In[2]:  dz_sc
       // Parameters_In[3]:  num_sc
       // Parameters_In[4]:  dz_vac  (is not used in this function)
       // Parameters_In[5]:  num_vac (is not used in this function)
       // Parameters_In[6..6+num_sc]: phi-array for the integration, Values in eV
       // Parameters_in[6+num_sc+1..6+num_sc+num_vac]: vacuum-barrier, Values in eV
       // Parameters_in array is created in I_tunnel
       function D_CB(W:double; Parameters_In:TFunctionParameters):double;
       var Parameters_out: TFunctionParameters;
           z_max,dz_sc:double;
           num_sc,i:  integer;
           voltage: double;
           Intersec: integer;
       begin
         voltage   :=   Parameters_In[0];
         z_max     :=   Parameters_In[1];
         dz_sc     :=   Parameters_In[2];
         num_Sc    :=   trunc(Parameters_In[3]);
         setlength(Parameters_out,num_sc+3);
         Parameters_out[0]:= W;       // W is in eV
         Parameters_out[1]:= dz_sc;   // dz in nm
         Parameters_out[2]:= num_sc;  // num is the size of the following array, which contains the band bedning for the integration
         for i := 0 to num_Sc-1 do
           Parameters_out[i+3]:=Parameters_in[i+6];

         //find intersection point of W and E_CB(z)
         Intersec:=0;
         while (Intersec<num_sc) and (SC.E_g+Parameters_in[Intersec+6]-W>0) do
           inc(Intersec);

         result:=exp(-2*sqrt(2*SC.m_cb*me*eVToJ/(h_bar*h_bar))*integrate(Sqrt_E_cb_minus_W,Parameters_out,0,Intersec*dz_sc,Integration_Accuracy)*1e-9);
       end;

// 2. Tunneling out of valence band thru space charge region (downward band bending)

       // 2a. Helping function for the integration of the barrier-potential in space charge region of the valence band
       // given by conduction band edge as is case 1. above
       // Note: Feenstras barriers is different and makes no sense
       // z: distance from SC surface
       // Parameters[0]: W
       // Parameters[1]: dz
       // Parameters[2]: num //Number of following array entries
       // Parameters[3..3+num]: phi-array for the integration as prepared in D_VB, Values in eV
       function Sqrt_W_minus_E_vb(z:double; Parameters: TFunctionParameters):double;
       var Re:double;
           num,i:integer;
           dz,iz:double;
           W:double;
           phi_x:double;
           E:double;
       begin
         W:=Parameters[0];
         dz:=Parameters[1];
         num:=trunc(Parameters[2]);

         //Linear interpolation of the potential for smooth integration
         iz:=z/dz;
         i:=trunc(iz);
         if i<num-1 then
           begin
             iz:=iz-i;
             phi_x:= Parameters[i+3]*(1-iz)+Parameters[(i+1)+3]*iz;
           end
         else
           phi_x:=Parameters[num-1+3];

         E:=SC.E_G+Phi_x;
         // if electrons fill the conduction band up to E_F then barrier for
         // electrons tunneling thru SCR increases to E_F since no free DOS
         // availaible at E_C:
         if SC.Inv_C then E:=max(E,SC.E_F);


         Re := (E-W); //(E_g+phi_x-W);    // E_CB + Phi - W
         if (Re>0) and (W>Phi_x) then
           result:=sqrt( Re )
         else
           result:=0;

       end;

       // 2b. calculates the transmission coefficient for tunneling thru the space charge region of the
       // valence band. see Eq. 4, 5 and Section E in [2]
       // W:        Energy in eV
       // Parameters_in[0]:  Voltage in V
       // Parameters_in[1]:  z_max_sc
       // Parameters_in[2]:  dz_sc
       // Parameters_in[3]:  num_sc
       // Parameters_in[4]:  dz_vac  (is not used in this function)
       // Parameters_in[5]:  num_vac (is not used in this function)
       // Parameters_in[6..6+num_sc]: phi-array for the integration, Values in eV
       // Parameters_in[6+num_sc+1..6+num_sc+num_vac]: vacuum-barrier, Values in eV
       // Parameters_in array is created in I_tunnel
       function D_VB(W:double; Parameters_In:TFunctionParameters):double;
       var Parameters_out: TFunctionParameters;
           z_max,dz_sc:double;
           D_vac:double; //transmission in vacuum
           D_SC: double; //transmission in space-charge region
           num_sc,i:  integer;
           voltage: double;
           Intersec: integer;
           phi_x,E:double;
       begin
         voltage   :=   Parameters_In[0];
         z_max     :=   Parameters_In[1];
         dz_sc     :=   Parameters_In[2];
         num_Sc    :=   trunc(Parameters_In[3]);

         setlength(Parameters_out,num_sc+3);
         Parameters_out[0]:= W;        // W is in eV
         Parameters_out[1]:= dz_sc;    // dz in nm
         Parameters_out[2]:= num_sc;   // num is the size of the following array, which contains the band bending for the integration
         for i := 0 to num_sc-1 do
           Parameters_out[i+3]:=Parameters_in[i+6];

         //find intersection point of W and E_VB(z)
         Intersec:=-1;
         repeat
           inc(Intersec);
           phi_x:=Parameters_in[Intersec+6];
           E:=SC.E_G+Phi_x;
           if SC.Inv_C then E:=max(E,SC.E_F);
         until not ((Intersec<num_sc) and (W>Phi_x));// and (SC.E_G+Phi_x-W>0));

         result:=exp(-2*sqrt(2*SC.m_vb*me*eVToJ/(h_bar*h_bar))*integrate(Sqrt_W_minus_E_vb,Parameters_out,0,Intersec*dz_sc,Integration_Accuracy)*1e-9);
       end;


// 3. Transmission coefficient for tunneling thru the vacuum barrier:
       // see Eq. 4, 5 and Section E in [2]
       // W:        Energy in eV
       // Parameters_in[0]:  Voltage in V
       // Parameters_in[1]:  z_max_sc
       // Parameters_in[2]:  dz_sc
       // Parameters_in[3]:  num_sc
       // Parameters_in[4]:  dz_vac
       // Parameters_in[5]:  num_vac
       function D_vac(W:double; Parameters_In:TFunctionParameters):double;
       var num_sc,num_vac:integer;
           dz_vac:double;
           voltage: double;
           sum:double;
           z:double;
           i: Integer;
           Vtrap:double;
       begin

         voltage:=   Parameters_In[0];
         num_sc  :=  trunc(Parameters_In[3]);
         num_vac :=  trunc(Parameters_In[5]);
         dz_vac  :=  Parameters_In[4];
         sum:=0;
         z:=0;
         for i := 0 to num_vac-1 do
           begin
             z:=i*dz_vac;
             Vtrap:= Parameters_In[6+num_sc+i];
             if Vtrap-W>0 then
               sum:=sum+sqrt(Vtrap-W);
           end;
         sum:=sum*dz_vac;
         result:=exp(-2*sqrt(2*me*eVToJ/(h_bar*h_bar))*sum*1e-9);

       end;


// 4. Multiplication of vacuum and space charge transmission coefficients:

       // tunneling into conduction band
       // Parameters_in is the same as in D_cb / D_vac
       function D_vac_times_D_cb(W:double;Parameters_In:TFunctionParameters):double;
       begin
         result:=D_vac(W,Parameters_In)*D_cb(W,Parameters_in);
       end;

       //tunneling out of valence band
       // Parameters_in is the same as in D_vb / D_vac
       function D_vac_times_D_vb(W:double;Parameters_In:TFunctionParameters):double;
       begin
         result:=D_vac(W,Parameters_In)*D_vb(W,Parameters_in);
       end;


//***************************************************************************************************************
//*              CALCULATION OF THE TUNNELING CURRENT COMPONENTS FOR THE CONDUCTION BAND                        *
//***************************************************************************************************************

       //Helping function for the integration of theta function times the integral of the transmission coefficient
       // E:        Energy in eV
       // Parameters[0]:  Voltage in V
       // Parameters[1]:  z_max_sc
       // Parameters[2]:  dz_sc
       // Parameters[3]:  num_sc
       // Parameters[4]:  dz_vac  (is not used in this function)
       // Parameters[5]:  num_vac (is not used in this function)
       // Parameters[6..6+num_sc]: phi-array for the integration, Values in eV
       // Parameters[6+num_sc+1..6+num_sc+num_vac]: vacuum-barrier, Values in eV
       // Parameters array is created in I_tunnel
       function Theta_times_int_D_CB_D_Vac_dW(E:double;Parameters:TFunctionParameters):double;
       var W_From, W_To, phi: double;
       begin
        //phi:=Parameters[6];
        //integration starts at E_CB (non degenerate doping) or at E_F (degenerate doping, i.e., E_F > E_CB)
        //Note: W is the normal energy only ( W= E - h^2*k_parallel^2/(2*m) )
        W_From:=E*(1-SC.m_cb)+SC.m_cb*E_CB_min;
        W_To  :=E;

        // Implementation of theta function times integration of the transmission coefficient,
        // see Eq. 2 in [2];
        if (E-E_CB_min) > 0 then
           result := integrate(D_vac_times_D_cb,Parameters,W_From,W_To,Integration_Accuracy)
        else
           result := 0;
       end;


       //***********************************************************************
       // 1. tunneling directly into and out of conduction band states,
       //     including space charge regions.
       //***********************************************************************
       function I_CB(Voltage:double;Parameters:TFunctionParameters):double;
       var phi:double;
       begin
         phi:=Parameters[6];
         result:=me*e*eVToJ*eVToJ/(2*Pi*Pi*power(h_bar,3))*integrate(Theta_times_int_D_CB_D_Vac_dW,Parameters,E_FQ_C,SC.E_F+Voltage,Integration_Accuracy);
         //result is in A/m^2 and has to be multiplied by the tunneling area
         result:=result*Area;
         output(format('I_cb = %-2.5e A (directly in/out of CB, including SCR)',[result]));
       end;




//***************************************************************************************************************
//*                 CALCULATION OF THE TUNNELING CURRENT COMPONENTS FOR THE VALENCE BAND                        *
//***************************************************************************************************************

//***************************************************************************************************************
// Tunneling out of the valence band thru space charge region and vacuum into the tip
// Voltage < 0 ; downward band bending (Phi < 0)
//***************************************************************************************************************

       // Helping function for the integration of theta function times the integral of the transmission coefficient
       // E:        Energy in eV
       // Parameters[0]:  Voltage in V
       // Parameters[1]:  z_max_sc
       // Parameters[2]:  dz_sc
       // Parameters[3]:  num_sc
       // Parameters[4]:  dz_vac  (is not used in this function)
       // Parameters[5]:  num_vac (is not used in this function)
       // Parameters[6..6+num_sc]: phi-array for the integration, Values in eV
       // Parameters[6+num_sc+1..6+num_sc+num_vac]: vacuum-barrier, Values in eV
       // Parameters array is created in I_tunnel
       function Theta_times_int_D_VB_D_vac_dW(E:double;Parameters:TFunctionParameters):double;
       var W_From, W_To, phi: double;
       begin
        //phi:=Parameters[6];
        //Note: W is the normal energy only ( W= E + h^2*k_parallel^2/(2*m) )
        //      In contrast to Feenstra: Since the energy of the valence band is defined such,
        //      that smaller values correspond to HIGHER energies, one has to add the parallel
        //      component of the energy h^2*k_parallel^2/(2*m).
        //      Just keep in mind: W is the normal (perpendicular to the surface) component
        //                         of the total energy E. Hence W can only reach values between
        //                         E_VB (whole energy is stored in parallel component)
        //                         and E (whole energy is stored in perpendicular component).
        //                         W = E_VB ... E.
        // Note: If you substitude k_parallel_max^2 = 2*m_eff*m/h_bar^2 * (E_VB - E) in
        //       W= E + h_bar^2*k_parallel^2/(2*m), one ends up in:
        //       W_max =  E + m_eff *(E_VB - E)
        //       which is equal to:
        //       W_max = E*(1-m_eff)+E_VB
        //       This is the boundary of the integration and is in priciple equal to the
        //       boundary of the integration for the conduction band.
        W_From:=E*(1-SC.m_vb)+SC.m_vb*E_VB_max;
        W_To  :=E;

        // Implementation of theta function times integration of the transmission coefficient,
        // see Eq. 2 in [2];
        if (E_VB_max-E)>0 then
            result := -integrate(D_vac_times_D_vb,Parameters,W_From,W_To,Integration_Accuracy)
         else
           result := 0;

       end;


       //**********************************************************************
       // 1a. Tunneling directly from valence band without space charge region
       // energy range for integration: min(E_VB + Phi_Surface,E_F)...E_F + Voltage
       //**********************************************************************

       function I_VB(Voltage:double;Parameters:TFunctionParameters):double;
       var phi:double;
       begin
         phi:=Parameters[6];
         result:=me*e*eVToJ*eVToJ/(2*Pi*Pi*power(h_bar,3))*integrate(Theta_times_int_D_VB_D_vac_dW,Parameters,E_FQ_V,SC.E_F+Voltage,Integration_Accuracy);
         //result is in A/m^2 and has to be multiplied by the tunneling area
         result:=result*Area;
         output(format('I_vb = %-2.5e A (directly out of/in VB, including SCR)',[result]));
       end;


//***************************************************************************************************************

       function InterpolateVtrap(Phi_central:array of TSpaceChargePoint; pz_array:array of double; z:double):double;
       var i,i_min,i_max:integer;
           n:integer;
           z_min,z_max:double;
           z0:double;
       begin
          n:=length(pz_array);

          i:=0;
          while (i<n-1) and (Phi_central[i].Material=mMetal) do
           inc(i);
          z0:=pz_array[i-1];

          if z0+z>0 then
          begin
            result:=0;
            exit;
          end;

          if z0+z<-d then
          begin
            result:=0;
            exit;
          end;

          i:=0;
          while (i<n-1) and (z>pz_array[i]-z0) do
           inc(i);

          if z0+z=pz_array[i] then
           begin
             result:=Phi_central[i].Phi;
           end
          else if i>0 then
          begin
             z_min:=pz_array[i-1];
             z_max:=pz_array[i];
             result:=Phi_central[i-1].Phi+((z0+z)-z_min)/(z_max-z_min)*(Phi_central[i].Phi-Phi_central[i-1].Phi);
          end
          else
            result:=Phi_central[0].Phi;
       end;



       function InterpolatePhi(Phi_central:array of TSpaceChargePoint; pz_array:array of double; z:double):double;
       var i,i_min,i_max:integer;
           n:integer;
           z_min,z_max:double;
       begin
          n:=length(pz_array);
          if z>pz_array[n-1] then
            begin
              result:=0;
              exit;
            end;
          i:=0;
          while (i<n-1) and (z>pz_array[i]) do
           inc(i);

          if z=pz_array[i] then
           begin
             result:=Phi_central[i].Phi;
           end
          else if i>0 then
          begin
             z_min:=pz_array[i-1];
             z_max:=pz_array[i];
             result:=Phi_central[i-1].Phi+(z-z_min)/(z_max-z_min)*(Phi_central[i].Phi-Phi_central[i-1].Phi);
          end
          else
            result:=Phi_central[0].Phi;
       end;




       //Calculates the tunneling current for the applied voltage.
       // Parameters:
       // Voltage: Voltage between tip and sample
       // BB:   an array containing all information of the band bending for the selected voltage:
       //       Structure is the same as the structure of Parameters:
       //          BB[0]:  Voltage in V
       //          BB[1]:  phi_s in eV
       //          BB[2]:  z_max in m
       //          BB[3]:  dz in m
       //          BB[4]:  num (the size of the following array, which contains the band bending for the integration in space charge region)
       //          BB[5..5+num]: band bending (phi_z) as a function of the position z (phi_z is sometimes denoted by phi_x, sorry)
       // result: sum over the single tunneling current contributions
       function I_tunnel(Voltage:double;  Phi_central:array of TSpaceChargePoint;
                         pz_array:array of double; var ivb,icb:double;
                         Semiconductor:TSemiconductor):double;
       var Parameters:TFunctionParameters;
           n,i_sc,i_vac:integer;
           sc_surf:integer;
           dz_sc:double;
           z,dz_vac:double;
           z_max:double;
           i:integer;
           num_sc:integer;
           num_vac:integer;
           lambda,Pot_mirror:double;
       begin
         Area:=100e-20; //in 1/m²
         SC:=Semiconductor;

         //Linear interpolation of band bending in semiconductor
         n:=length(Phi_Central);
         i_vac:=0;
         i_sc:=0;


         E_VB_Max:=0;
         E_CB_Min:=0;

         for i := 0 to n-1 do
           begin
             if Phi_Central[i].Material=mVac then inc(i_vac);
             if Phi_Central[i].Material=mSC then
             begin
               inc(i_Sc);
               E_CB_Min:=min(E_CB_Min,Phi_Central[i].Phi);
               E_VB_Max:=max(E_VB_Max,Phi_Central[i].Phi);
             end;
             if (i>0) and (Phi_Central[i-1].Material=mVac) and (Phi_Central[i].Material=mSC) then
               sc_surf:=i;
           end;

         E_CB_Min:=E_CB_Min+SC.E_g;

         num_sc:=0;
         z_max:=0;
         dz_sc:=stepwidth_z;

         while (z_max<pz_array[length(pz_array)-1])  do
           begin
             inc(num_sc);
             setlength(Parameters,num_sc+6);
             Parameters[num_sc+6-1]:=InterpolatePhi(Phi_central,pz_array,z_max);
             z_max:=z_max+dz_sc;
           end;
         z_max:=z_max-dz_sc;

         z_max:=min(z_max,z_stop);
         //Linear interpolation of potential in Vacuum

         num_vac:=1000;
         dz_vac:=d/(num_vac-1);
         setlength(Parameters,num_sc+num_vac+6);
         lambda:= e*e*ln(2)/(8*Pi*epsi_0*d*1e-9);

         Parameters[num_sc+6]:= -1e10;
         Parameters[num_sc+6+num_vac-1]:= -1e10;

         for i := 1 to num_vac-2 do
           begin
             z:=i*dz_vac;
             Pot_mirror:= (1.15/2*lambda*d*d)/(z*(d-z))/eVToJ;
             Parameters[num_sc+6+i]:=InterpolateVTrap(Phi_central,pz_array,z)+SC.E_g+SC.Chi-Pot_mirror;
           end;

         //       Parameters[0]:  Voltage in V
         //       Parameters[1]:  z_max_sc
         //       Parameters[2]:  dz_sc
         //       Parameters[3]:  num_sc
         //       Parameters[4]:  dz_vac
         //       Parameters[5]:  num_vac


         Parameters[0]:=Voltage;
         Parameters[1]:=z_max;
         Parameters[2]:=dz_sc;
         Parameters[3]:=num_sc;
         Parameters[4]:=dz_vac;
         Parameters[5]:=num_vac;

         if I_V_output then
         begin
           output(format('\nCalculation of the tunneling current (BIAS = %.5e V):\n',[voltage]));
           ivb:=I_VB(Voltage,Parameters);
           icb:=I_CB(Voltage,Parameters);
           result := ivb + icb;
         end
         else
           begin
            ivb:=0;
            icb:=0;
            result := 0;
           end;





       end;



end.
