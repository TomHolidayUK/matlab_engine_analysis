clear
clc
 
%MECH3496 - Thermofluids III - Aerospace Propulsion Assignment
%Tom Holiday
%201015193

%%%%This programme returns a full set of working gas properties for a General
%Electric CFM56-5C non-mixing turbofan engine by implementing a thermodynamic cycle analysis.
%The input parameters can be changed to model any non-mixing turbofan
%engine. The model goes on to calculate thrust and TSFC and compare the results over a rnage of Turbine Entry Temperatures (TET's) This model makes assumptions which are listed below

%%%Analysis Model Assumptions
%The intake is assumed to have no temperature loss and to be adiabatic
%(Tt1 = Tt0) 
%%Jet Pipe assumed to be Adiabatic (Tt7 = Tt6)
%Core nozzle assumed to be Adiabatic (Tt8 = Tt7)
%Bypass duct assumed to be Adiabatic (Tt9 = Tt2)

%%%Stages
%(0-1) Intake Stage
%(1-2) LP compressor stage
%(2-3) HP compressor stage 
%(3-4) Combustor stage 
%(4-5) HP turbine stage 
%(5-6) LP turbine stage 
%(6-7) Jet pipe stage 
%(7-8) Core (hot jet) nozzle stage 
%(2-9) Bypass (cold jet) nozzle stage 

%%%Data Notation
%Pn = Static Pressure at Point n
%Ptn = Total Pressure at Point n
%Tn = Static Temperature at Point n
%Ttn = Total Temperature at Point n
%Ttn = Total Isentropic Temperature at Point n
%Vn = Air velocity at Point n
 
% Bypass Ratio (BPR)
BPR = 7.41;
% Normalised Spool Speed (N)
N = 1.19;

%Input Parameters

P0 = 19416.7; %Ambient Pressure (from ISA)
T0 = 216.65; %Ambient Temperature (from ISA)
M0 = 0.8; % Intake Mach Number
ycold = 1.4; % Cold heat capacity ratio 
Rcold = 287.0; % Cold gas constant 
PRi = 0.96; % Intake pressure recovery 
PRbypass = 0.9; % Bypass nozzle pressure recovery

% (0-1) Intake Stage
 
V0 = M0 * sqrt(ycold * Rcold * T0);

Pt0 = P0 * (1 + ((ycold - 1) / 2) * M0^2)^(ycold / (ycold -1 ));

Tt0 = T0 * (Pt0 / P0)^((ycold - 1)/ycold);

Pt1 = PRi*Pt0;

Tt1 = Tt0;

%Calculating PRlpc for when bypass nozzle is just choked

PRcritbypass = ((ycold + 1) / 2)^(ycold/(ycold-1));

PRlpc = (PRcritbypass * P0) / (PRbypass * Pt1); % LP compressor pressure ratio (as calculated in Q2a)

% Turbine entry temperature (TET) (n), ranging from n1 to n3 by an
% increment of n2. The user can change these as they wish.
 
n = 1;
n1 = 1050;
n2 = 1;
n3 = 2500;

for TET = n1:n2:n3;
 
%Further Input Parameters

MFRair = N*408; % Intake air mass flow rate 
IElpc = 0.93; % LP compressor isentropic efficiency  
IEhpc = 0.91; % HP compressor isentropic efficiency
PRcc = 0.98; % Combustor pressure ratio 
Ecc = 0.99; % Combustor efficiency 
LHV = 43.1*10^6; % Fuel lower heating value 
IEhpt = 0.94; % HP turbine isentropic efficiency 
MEhpt = 0.99; % HP spool mechanical efficiency 
IElpt = 0.95; % LP turbine isentropic efficiency 
MElpt = 0.99; % LP spool mechanical efficiency 
PRjp = 1.00; % Jet pipe pressure recovery 
PRcn = 1.00; % Core nozzle pressure recovery 
yhot = 1.333; % Hot heat capacity ratio 
Cpcold = 1005; % Cold specific heat at constant pressure 
Cphot = 1150; % Hot specific heat at constant pressure
Rhot = 287.3; % Hot gas constant 

%%%%(1-2) LP Compressor Stage
 
Pt2 = PRlpc * Pt1;

Tt2i = Tt1 * (Pt2/Pt1)^((ycold-1)/ycold);

Tt2 = ((Tt2i-Tt1) / IElpc) + Tt1;
 
%%%%%(2-3) HP compressor stage 

PRhpc = 26.7 * (N^1.22); % HP compressor pressure ratio
 
Pt3 = PRhpc * Pt2;

Tt3i = Tt2 * (Pt3 / Pt2)^((ycold-1)/ycold);

Tt3 = ((Tt3i-Tt2) / IEhpc) + Tt2;
 
%(3-4) Combustor stage 
 
Tt4 = TET;

Pt4 = PRcc * Pt3;

AFR = ((LHV * Ecc)/(Cphot * (Tt4 - Tt3)))-1;
 
%(4-5) HP turbine stage 
 
Tt5 = Tt4 - (AFR * Cpcold * (Tt3 - Tt2))/((AFR + 1) * Cphot * MEhpt);

Tt5i = Tt4 - (Tt4 - Tt5) / IEhpt;

Pt5 = Pt4 * (Tt5i / Tt4)^(yhot/(yhot-1));
 
%(5-6) LP turbine stage 
 
MFRcore = MFRair / (1+BPR);

Tt6 = Tt5 - ((1 + BPR) * AFR * Cpcold * (Tt2 - Tt1))/((AFR + 1) * Cphot * MElpt);

Tt6i = Tt5 - (Tt5 - Tt6) / IElpt;

Pt6 = Pt5 * (Tt6i / Tt5)^(yhot/(yhot-1));
 
%(6-7) Jet pipe stage 

Pt7 = Pt6 * PRjp;

Tt7 = Tt6; %Jet Pipe assumed to be Adiabatic
 
%(7-8) Core (hot jet) nozzle stage 

%No total pressure or total temperature loss at core nozzle as there is no afterburner.
 
Pt8 = Pt7 * PRcn;

Tt8 = Tt7; %Core nozzle assumed to be Adiabatic

PRcritcore = ((yhot + 1) / 2)^(yhot/(yhot-1)); %Critical Pressure of Core

PRidealcore = Pt8 / P0; %Ideal Pressure Ratio of Core

TRcritcore = (yhot + 1) / 2; %Critical Temperature Ratio of Core

P8 = Pt8 / PRcritcore;

T8 = Tt8 / TRcritcore; 

if PRidealcore >= PRcritcore
    M8 = 1;
    M8_Choked_Condition = "Choked";
else
    M8 = (sqrt(2 * Cphot * (Tt8 - T8))) / (sqrt(yhot * Rhot * T8));
    M8_Choked_Condition = "Not Choked";
end

V8 = M8 * sqrt(yhot * Rhot * T8);

A8 = ((1 + (1 / AFR)) * ((Rhot * T8) / (P8 * V8))) * MFRcore;

Fncore = ((1 + (1 / AFR)) * V8 - V0 + ((A8 / MFRcore) * (P8 - P0))) * MFRcore;
 
%(2-9) Bypass (cold jet) nozzle stage 
 
Tt9 = Tt2; %Bypass duct assumed to be Adiabatic

Pt9 = PRbypass * Pt2;

PRcritbypass = ((ycold + 1) / 2)^(ycold/(ycold-1));

PRidealbypass = Pt9 / P0;

TRcritbypass = (ycold + 1) / 2;

P9 = Pt9 / PRcritbypass;

T9 = Tt9 / TRcritbypass;

if PRidealbypass >= PRcritbypass
    M9 = 1;
    M9_Choked_Condition = "Choked";
elseif abs(PRidealbypass - PRcritbypass) < 0.001 %This is essentially 'PRidealbypass == PRcritbypass' but puts a tolerance of 0.001 in so they dont have to be EXACTLY equal 
    M9 = 1;
    M9_Choked_Condition = "Just Choked";
else
    M9 = (sqrt(2 * Cpcold * (Tt9 - T9))) / (sqrt(ycold * Rcold * T9));
    M9_Choked_Condition = "Not Choked";
end

V9 = M9 * sqrt(ycold * Rcold * T9);
 
A9 = (BPR * (Rhot * T8)/(P8 * V8)) * MFRcore;

Fnbypass = (BPR * (V9 - V0) + (A9 / MFRcore) * (P9 - P0)) * MFRcore;
 
%%%Thrust Calculations%%%

Fntotal = Fncore + Fnbypass;
 
Fntotalvector(n) = Fntotal;
 
MFRbypass = MFRcore * BPR;

MFRfuel = MFRcore / AFR;
 
TSFC = MFRfuel / Fntotal; % Total specific fuel consumption
 
TSFCvector(n) = TSFC;

%The code below checks to make sure the thrust generated is positive for the complete range of TET values, if
%not (as a results of TET beign too low) it tells the user to increase n1 -
%the lower limit of TET.

Fn_Thrust_Condition = Fntotalvector(1);

if Fn_Thrust_Condition < 0
   Thrust_Condition = "Negative Thrust - Invalid Results - Run again with higher minimum TET (n1)";
else
   Thrust_Condition = "Positive Thrust - Valid Results";
end



%-------------------------------------------------------------------------



%%%Sensitivity Study 1

PRhpc_SS1 = PRhpc * 1.1; %When PRhpc increases by 10%

%All parameters downstream of the HP Comressor now different (marked _SS1)

Pt3_SS1 = PRhpc_SS1 * Pt2; 

Tt3i_SS1 = Tt2 * (Pt3_SS1 / Pt2)^((ycold-1)/ycold);

Tt3_SS1 = ((Tt3i_SS1 - Tt2) / IEhpc) + Tt2;
 
%(3-4) Combustor stage 
 
Tt4 = TET;

Pt4_SS1 = PRcc * Pt3_SS1;

AFR_SS1 = ((LHV * Ecc)/(Cphot * (Tt4 - Tt3_SS1))) - 1;
 
%(4-5) HP turbine stage 
 
Tt5_SS1 = Tt4 - (AFR_SS1 * Cpcold * (Tt3_SS1 - Tt2))/((AFR_SS1 + 1) * Cphot * MEhpt);

Tt5i_SS1 = Tt4 - (Tt4 - Tt5_SS1) / IEhpt;

Pt5_SS1 = Pt4_SS1 * (Tt5i_SS1 / Tt4)^(yhot/(yhot-1));
 
%(5-6) LP turbine stage 
 
MFRcore = MFRair / (1 + BPR);

Tt6_SS1 = Tt5_SS1 - ((1 + BPR) * AFR_SS1 * Cpcold * (Tt2 - Tt1))/((AFR_SS1 + 1) * Cphot * MElpt);

Tt6i_SS1 = Tt5_SS1 - (Tt5_SS1 - Tt6_SS1) / IElpt;

Pt6_SS1 = Pt5_SS1 * (Tt6i_SS1 / Tt5_SS1)^(yhot/(yhot-1));
 
%(6-7) Jet pipe stage 

Pt7_SS1 = Pt6_SS1 * PRjp;

Tt7_SS1 = Tt6_SS1; %Jet Pipe assumed to be Adiabatic
 
%(7-8) Core (hot jet) nozzle stage 
 
Pt8_SS1 = Pt7_SS1 * PRcn;

Tt8_SS1 = Tt7_SS1; %Core Nozzle assumed to be Adiabatic

PRcritcore = ((yhot + 1) / 2)^(yhot/(yhot-1)); %Critical Pressure of Core

PRidealcore_SS1 = Pt8_SS1 / P0; %Ideal Pressure Ratio of Core

TRcritcore = (yhot + 1) / 2; %Critical Temperature Ratio of Core

P8_SS1 = Pt8_SS1 / PRcritcore;

T8_SS1 = Tt8_SS1 / TRcritcore;

if PRidealcore_SS1 >= PRcritcore
    M8_SS1 = 1;
    M8_SS1_Choked_Condition = "Choked";
else
    M8_SS1 = (sqrt(2 * Cphot * (Tt8_SS1 - T8_SS1))) / (sqrt(yhot * Rhot * T8_SS1));
    M8_SS1_Choked_Condition = "Not Choked";
end

V8_SS1 = M8_SS1 * sqrt(yhot * Rhot * T8_SS1);

A8_SS1 = ((1 + (1 / AFR_SS1)) * ((Rhot * T8_SS1) / (P8_SS1 * V8_SS1))) * MFRcore;

Fncore_SS1 = ((1 + (1 / AFR_SS1)) * V8_SS1 - V0 + ((A8_SS1 / MFRcore) * (P8_SS1 - P0))) * MFRcore;
 
%(2-9) Bypass (cold jet) nozzle stage 
 
Tt9 = Tt2; %Bypass Duct assumed to be Adiabatic

Pt9 = PRbypass * Pt2;

PRcritbypass = ((ycold + 1) / 2)^(ycold/(ycold-1));

PRidealbypass = Pt9 / P0;

TRcritbypass = (ycold + 1) / 2;

P9 = Pt9 / PRcritbypass;

T9 = Tt9 / TRcritbypass;

if PRidealbypass >= PRcritbypass
    M9 = 1;
    M9_SS1_Choked_Condition = "Choked";
elseif abs(PRidealbypass - PRcritbypass) < 0.001 %This is essentially 'PRidealbypass == PRcritbypass' but puts a tolerance of 0.001 in so they dont have to be EXACTLY equal
    M9 = 1;
    M9_SS1_Choked_Condition = "Just Choked";
else
    M9 = (sqrt(2 * Cpcold * (Tt9 - T9))) / (sqrt(ycold * Rcold * T9));
    M9_SS1_Choked_Condition = "Not Choked";
end

V9 = M9 * sqrt(ycold * Rcold * T9);
 
A9_SS1 = (BPR * (Rhot * T8_SS1)/(P8_SS1 * V8_SS1)) * MFRcore;

Fnbypass_SS1 = (BPR * (V9 - V0) + (A9_SS1 / MFRcore) * (P9 - P0)) * MFRcore;
 
%%%Thrust Calculations%%%

Fntotal_SS1 = Fncore_SS1 + Fnbypass_SS1;
 
Fntotalvector_SS1(n) = Fntotal_SS1;
 
MFRbypass = MFRcore * BPR;

MFRfuel_SS1 = MFRcore / AFR_SS1;
 
TSFC_SS1 = MFRfuel_SS1 / Fntotal_SS1; % Total specific fuel consumption
 
TSFCvector_SS1(n) = TSFC_SS1;

%The code below checks to make sure the thrust generated is positive (for SS1) for the complete range of TET values, if
%not (as a results of TET beign too low) it tells the user to increase n1 -
%the lower limit of TET.

Fn_Thrust_Condition_SS1 = Fntotalvector_SS1(1);

if Fn_Thrust_Condition_SS1 < 0
   Thrust_Condition_SS1 = "Negative Thrust - Invalid Results - Run again with higher minimum TET (n1)";
else
   Thrust_Condition_SS1 = "Positive Thrust - Valid Results";
end


%-------------------------------------------------------------------------



%%%Sensitivity Study 2

PRhpc_SS2 = PRhpc * 0.9; %When PRhpc decreases by 10%

%All parameters downstream of the HP Comressor now different (marked _SS2)

Pt3_SS2 = PRhpc_SS2 * Pt2; 

Tt3i_SS2 = Tt2 * (Pt3_SS2 / Pt2)^((ycold-1)/ycold);

Tt3_SS2 = ((Tt3i_SS2 - Tt2) / IEhpc) + Tt2;
 
%(3-4) Combustor stage 
 
Tt4 = TET;

Pt4_SS2 = PRcc * Pt3_SS2;

AFR_SS2 = ((LHV * Ecc)/(Cphot * (Tt4 - Tt3_SS2))) - 1;
 
%(4-5) HP turbine stage 
 
Tt5_SS2 = Tt4 - (AFR_SS2 * Cpcold * (Tt3_SS2 - Tt2))/((AFR_SS2 + 1) * Cphot * MEhpt);

Tt5i_SS2 = Tt4 - (Tt4 - Tt5_SS2) / IEhpt;

Pt5_SS2 = Pt4_SS2 * (Tt5i_SS2 / Tt4)^(yhot/(yhot-1));
 
%(5-6) LP turbine stage 
 
MFRcore = MFRair / (1 + BPR);

Tt6_SS2 = Tt5_SS2 - ((1 + BPR) * AFR_SS2 * Cpcold * (Tt2 - Tt1))/((AFR_SS2 + 1) * Cphot * MElpt);

Tt6i_SS2 = Tt5_SS2 - (Tt5_SS2 - Tt6_SS2) / IElpt;

Pt6_SS2 = Pt5_SS2 * (Tt6i_SS2 / Tt5_SS2)^(yhot/(yhot-1));
 
%(6-7) Jet pipe stage 

Pt7_SS2 = Pt6_SS2 * PRjp;

Tt7_SS2 = Tt6_SS2; %Jet Pipe assumed to be Adiabatic
 
%(7-8) Core (hot jet) nozzle stage 
 
Pt8_SS2 = Pt7_SS2 * PRcn;

Tt8_SS2 = Tt7_SS2; %Core Nozzle assumed to be Adiabatic

PRcritcore = ((yhot + 1) / 2)^(yhot/(yhot-1)); %Critical Pressure of Core

PRidealcore_SS2 = Pt8_SS2 / P0; %Ideal Pressure Ratio of Core

TRcritcore = (yhot + 1) / 2; %Critical Temperature Ratio of Core

P8_SS2 = Pt8_SS2 / PRcritcore;

T8_SS2 = Tt8_SS2 / TRcritcore;

if PRidealcore_SS2 >= PRcritcore
    M8_SS2 = 1;
    M8_SS2_Choked_Condition = "Choked";
else
    M8_SS2 = (sqrt(2 * Cphot * (Tt8_SS2 - T8_SS2))) / (sqrt(yhot * Rhot * T8_SS2));
    M8_SS2_Choked_Condition = "Not Choked";
end

V8_SS2 = M8_SS2 * sqrt(yhot * Rhot * T8_SS2);

A8_SS2 = ((1 + (1 / AFR_SS2)) * ((Rhot * T8_SS2) / (P8_SS2 * V8_SS2))) * MFRcore;

Fncore_SS2 = ((1 + (1 / AFR_SS2)) * V8_SS2 - V0 + ((A8_SS2 / MFRcore) * (P8_SS2 - P0))) * MFRcore;
 
%(2-9) Bypass (cold jet) nozzle stage 
 
Tt9 = Tt2; %Bypass Duct assumed to be Adiabatic

Pt9 = PRbypass * Pt2;

PRcritbypass = ((ycold + 1) / 2)^(ycold/(ycold-1));

PRidealbypass = Pt9 / P0;

TRcritbypass = (ycold + 1) / 2;

P9 = Pt9 / PRcritbypass;

T9 = Tt9 / TRcritbypass;

if PRidealbypass >= PRcritbypass
    M9 = 1;
    M9_SS2_Choked_Condition = "Choked";
elseif abs(PRidealbypass - PRcritbypass) < 0.001 %This is essentially 'PRidealbypass == PRcritbypass' but puts a tolerance of 0.001 in so they dont have to be EXACTLY equal
    M9 = 1;
    M9_SS2_Choked_Condition = "Just Choked";
else
    M9 = (sqrt(2 * Cpcold * (Tt9 - T9))) / (sqrt(ycold * Rcold * T9));
    M9_SS2_Choked_Condition = "Not Choked";
end

V9 = M9 * sqrt(ycold * Rcold * T9);
 
A9_SS2 = (BPR * (Rhot * T8_SS2)/(P8_SS2 * V8_SS2)) * MFRcore;

Fnbypass_SS2 = (BPR * (V9 - V0) + (A9_SS2 / MFRcore) * (P9 - P0)) * MFRcore;
 
%%%Thrust Calculations%%%

Fntotal_SS2 = Fncore_SS2 + Fnbypass_SS2;
 
Fntotalvector_SS2(n) = Fntotal_SS2;
 
MFRbypass = MFRcore * BPR;

MFRfuel_SS2 = MFRcore / AFR_SS2;
 
TSFC_SS2 = MFRfuel_SS2 / Fntotal_SS2; % Total specific fuel consumption
 
TSFCvector_SS2(n) = TSFC_SS2;

%The code below checks to make sure the thrust generated is positive (for SS2) for the complete range of TET values, if
%not (as a results of TET beign too low) it tells the user to increase n1 -
%the lower limit of TET.

Fn_Thrust_Condition_SS2 = Fntotalvector_SS2(1);

if Fn_Thrust_Condition_SS2 < 0
   Thrust_Condition_SS2 = "Negative Thrust - Invalid Results - Run again with higher minimum TET (n1)";
else
   Thrust_Condition_SS2 = "Positive Thrust - Valid Results";
end

%----------------------------------------------------------------------
 
n = n + 1;
 
end

TET = n1:n2:n3;

plot(TET,Fntotalvector)

hold on

plot(TET,Fntotalvector_SS1)

hold on

plot(TET,Fntotalvector_SS2)

xlabel('Turbine Entry Temperature, K')
ylabel('Thrust, N')
title('Thrust vs TET')
legend('Standard HP Pressure Ratio','10% Increase in HP Pressure Ratio','10% Decrease in HP Pressure Ratio')
 
figure(2)

plot(TET,TSFCvector)

hold on 

plot(TET,TSFCvector_SS1)

hold on 

plot(TET,TSFCvector_SS2)

xlabel('Turbine Entry Temperature, K')
ylabel('Thrust Specific Fuel Consumption, kg.s/N')
title('TSFC vs TET')
legend('Standard HP Pressure Ratio','10% Increase in HP Pressure Ratio','10% Decrease in HP Pressure Ratio')

%Finding minimum value of TSFC and corresponding value of TET

[TSFCmin,ni] = min(TSFCvector);

TETmin = 1025 + ni;

%Finding minimum value of TSFC and corresponding value of TET for when PRhpc increases by 10%

[TSFCmin_SS1,ni_SS1] = min(TSFCvector_SS1);

TETmin_SS1 = 1025 + ni_SS1;

%Finding minimum value of TSFC and corresponding value of TET for when PRhpc decreases by 10%

[TSFCmin_SS2,ni_SS2] = min(TSFCvector_SS2);

TETmin_SS2 = n1 + ni_SS2;

%end



