clear
clc
 
%MECH3496 - Thermofluids III - Aerospace Propulsion Assignment
%Tom Holiday
%201015193

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

% Turbine entry temperature (TET) of 1189 as this has been calculated to
% return a minimum value of TSFC.

TET = 1189;
 
%Engine Performance Data

MFRair = N*408; % Intake air mass flow rate 
PRi = 0.96; % Intake pressure recovery  
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
PRbypass = 0.9; % Bypass nozzle pressure recovery
ycold = 1.4; % Cold heat capacity ratio 
yhot = 1.333; % Hot heat capacity ratio 
Cpcold = 1005; % Cold specific heat at constant pressure 
Cphot = 1150; % Hot specific heat at constant pressure
Rcold = 287.0; % Cold gas constant 
Rhot = 287.3; % Hot gas constant 
P0 = 19416.7; %Ambient Pressure (from ISA)
T0 = 216.65; %Ambient Temperature (from ISA)
M0 = 0.8; % Intake Mach Number

% (0-1) Intake Stage
 
V0 = M0 * sqrt(ycold * Rcold * T0);

Pt0 = P0 * (1 + ((ycold - 1) / 2) * M0^2)^(ycold / (ycold -1 ));

Tt0 = T0 * (Pt0 / P0)^((ycold - 1)/ycold);

Pt1 = PRi*Pt0;

Tt1 = Tt0;

%Calculating PRlpc for when bypass nozzle is just choked

PRlpc = 1.4373; % LP compressor pressure ratio (calculated in Q2a)
 
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
 
Pt8 = Pt7 * PRcn;

Tt8 = Tt7; %Core nozzle assumed to be Adiabatic

PRcritcore = ((yhot + 1) / 2)^(yhot/(yhot-1)); %Critical Pressure of Core

PRidealcore = Pt8 / P0; %Ideal Pressure Ratio of Core

TRcritcore = (yhot + 1) / 2; %Critical Temperature Ratio of Core

P8 = Pt8 / PRcritcore;

T8 = Tt8 / TRcritcore; 

if PRidealcore >= PRcritcore
    M8 = 1;
    M8_Choked_Condition = "Yes";
else
    M8 = (sqrt(2 * Cphot * (Tt8 - T8))) / (sqrt(yhot * Rhot * T8));
    M8_Choked_Condition = "No";
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
elseif PRidealbypass == PRcritbypass
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
 
MFRbypass = MFRcore * BPR;

MFRfuel = MFRcore / AFR;
 
TSFC = MFRfuel / Fntotal % Total specific fuel consumption
 