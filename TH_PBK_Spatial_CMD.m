clear all
clc
tic
%--------- Time unit: second (s), concentration unit: pM, volume unit: L-------%


%% ------------------------- Default Parameters ----------------------------------------- %%
BW    = 75;
QTC   = 0.015;
QLC   = 0.25;
QRBC  = 1 - QTC - QLC;
VTC   = 0.0003;
VLC   = 0.0257;
VPC   = 0.0424;
VRBC  = 0.91 - VTC - VLC - VPC;
HCT   = 0.44;
VTBC  = 0.18;
VLBC  = 0.11;
VRBBC = 0.0236;
VT    = VTC*BW;
VRB   = VRBC*BW;
VL    = VLC*BW;
VPtot = VPC*BW; %Total plasma volume

param.QC   = 15*BW^0.74/3600*(1-HCT);  %Plasma cardiac output
param.QT   = QTC*param.QC;  %Plasma flow rate to Thyroid
param.QRB  = QRBC*param.QC; %Plasma flow rate to RB
param.QL   = QLC*param.QC;  %Plasma flow rate to Liver
param.VTB  = VT*VTBC*(1-HCT);   %Plasma volume in Thyroid
param.VRBB = VRB*VRBBC*(1-HCT); %Plasma volume in RB
param.VRBT = VRB*(1-VRBBC); %Tissue volume in RB
param.VLB  = VL*VLBC*(1-HCT);   %Plasma volume in Liver
param.VLT  = VL*(1-VLBC);   %Tissue volume in Liver
param.VB   = VPtot - param.VTB - param.VRBB - param.VLB;    %Plasma volume of Body Blood compartment

param.fuT4RBT = 0.1;
param.fuT3RBT = 0.01;
param.fuT4LT  = 0.1;
param.fuT3LT  = 0.01;

param.TBGtot = 3.515E5*VPtot;
param.TTRtot = 5.35E6*VPtot;
param.ALBtot = 6.45E8*VPtot;

param.kdT4TBG = 60;
param.kdT4TTR = 5000;
param.kdT4ALB = 1.33E6;
param.kdT3TBG = 1100;
param.kdT3TTR = 3.25E5;
param.kdT3ALB = 9.75E6;

param.k2 = 0.018;  
param.k4 = 0.0832; 
param.k6 = 1.3; 
param.k8 = 0.165; 
param.k10 = 0.69; 
param.k12 = 2.2;

param.k1 = param.k2/param.kdT4TBG;
param.k3 = param.k4/param.kdT4TTR;
param.k5 = param.k6/param.kdT4ALB;
param.k7 = param.k8/param.kdT3TBG;
param.k9 = param.k10/param.kdT3TTR;
param.k11 = param.k12/param.kdT3ALB;

param.k20 = 1.6; %pmol/s
ThyroidT4dailyprod = param.k20/1E12*776.87*1E6*24*3600; %ug/day. T4 MW=776.87.
param.k21 = 1;
param.k22 = 1.6/14; %pMol/s
ThyroidT3dailyprod = param.k22/1E12*651*1E6*24*3600; %ug/day. %T3 MW=651.
param.k23 = 0.15;

a1 = 0.25; 
a2 = 0.25;

param.k24 = 1.7084E-6/param.fuT4RBT*a2;
param.k25 = 4.243;
param.k26 = 1.44E-6/param.fuT4LT*a1;
param.k27 = 1.398;
param.k28 = 0.00141/param.fuT4RBT;
param.k29 = 8.9E-4/param.fuT3RBT;
param.k30 = 2.786E-4/param.fuT4LT;
param.k31 = 1.9688E-3/param.fuT3LT;
param.k32 = 1.7084E-6/param.fuT4RBT*(1-a2);
param.k33 = 6.3349E-6/param.fuT3RBT;
param.k34 = 1.44E-6/param.fuT4LT*(1-a1);
param.k35 = 3.4229E-5/param.fuT3LT; 

%Parameters for endocrine disrupting chemical X production and clearance in the Body Blood compartment(VB)
param.k37 = 0;
param.k36 = 0;

%Parameters for X binding to TTR
param.kdXTTR = 5000;
param.k39 = 0;%0.0832;
param.k38 = param.k39/param.kdXTTR;

%Parameters for X binding to TBG
param.kdXTBG = 60;
param.k43 = 0;%0.018;
param.k42 = param.k43/param.kdXTBG;

param.M = 20; %Number of tissue segments

default_param = struct();
fields = fieldnames(param);
for i = 1:numel(fields)
    default_param.(fields{i}) = param.(fields{i});
end

%% ------------------------- Initial Condition ------------------------------------------ %%
init.fT4B    = 15;
init.T4TBGB  = 0;
init.T4TTRB  = 0;
init.T4ALBB  = 0;
init.fT3B    = 5;
init.T3TBGB  = 0;
init.T3TTRB  = 0;
init.T3ALBB  = 0;
init.TBGB    = param.TBGtot/param.VB;
init.TTRB    = param.TTRtot/param.VB;
init.ALBB    = param.ALBtot/param.VB;

init.fT4T    = 15;
init.T4TBGT  = 0;
init.T4TTRT  = 0;
init.T4ALBT  = 0;
init.fT3T    = 5;
init.T3TBGT  = 0;
init.T3TTRT  = 0;
init.T3ALBT  = 0;
init.TBGT    = 0;
init.TTRT    = 0;
init.ALBT    = 0;

init.XB      = 0;
init.XTTRB   = 0;
init.XT      = 0;
init.XTTRT   = 0;
init.XTBGB   = 0;
init.XTBGT   = 0;

init.fT4RB(1:param.M,1)   = 15;
init.T4TBGRB(1:param.M,1) = 0;
init.T4TTRRB(1:param.M,1) = 0;
init.T4ALBRB(1:param.M,1) = 0;
init.fT3RB(1:param.M,1)   = 5;
init.T3TBGRB(1:param.M,1) = 0;
init.T3TTRRB(1:param.M,1) = 0;
init.T3ALBRB(1:param.M,1) = 0;
init.TBGRB(1:param.M,1)   = 0;
init.TTRRB(1:param.M,1)   = 0;
init.ALBRB(1:param.M,1)   = 0;
init.T4RBT(1:param.M,1)   = 0;
init.T3RBT(1:param.M,1)   = 0;
init.XRB(1:param.M,1)     = 0;
init.XTTRRB(1:param.M,1)  = 0;
init.XTBGRB(1:param.M,1)  = 0;

init.fT4L(1:param.M,1)    = 15;
init.T4TBGL(1:param.M,1)  = 0;
init.T4TTRL(1:param.M,1)  = 0;
init.T4ALBL(1:param.M,1)  = 0;
init.fT3L(1:param.M,1)    = 5;
init.T3TBGL(1:param.M,1)  = 0;
init.T3TTRL(1:param.M,1)  = 0;
init.T3ALBL(1:param.M,1)  = 0;
init.TBGL(1:param.M,1)    = 0;
init.TTRL(1:param.M,1)    = 0;
init.ALBL(1:param.M,1)    = 0;
init.T4LT(1:param.M,1)    = 0;
init.T3LT(1:param.M,1)    = 0;
init.XL(1:param.M,1)      = 0;
init.XTTRL(1:param.M,1)   = 0;
init.XTBGL(1:param.M,1)   = 0;

init_default = struct();
fields = fieldnames(init);
for i = 1:numel(fields)
    init_default.(fields{i}) = init.(fields{i});
end

%% ------------------------- Run To Steady-State (Table 1, S1-S4) ----------------------- %%
y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600];  %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode15s('TH_PBK_Spatial_ODE',tspan, y0, options, param); 

y0_default_steady_state = y(end, :); %Using the steady-state values from the "Run to Steady State" section above as the initial values

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.XB      = y(:,23);
model.XTTRB   = y(:,24);
model.XT      = y(:,25);
model.XTTRT   = y(:,26);
model.XTBGB   = y(:,27);
model.XTBGT   = y(:,28);

N = 29;  %Starting index of repeating variables
M = param.M;
model.fT4RB   = y(:,N:N+M-1);
model.T4TBGRB = y(:,N+M:N+2*M-1);
model.T4TTRRB = y(:,N+2*M:N+3*M-1);
model.T4ALBRB = y(:,N+3*M:N+4*M-1);
model.fT3RB   = y(:,N+4*M:N+5*M-1);
model.T3TBGRB = y(:,N+5*M:N+6*M-1);
model.T3TTRRB = y(:,N+6*M:N+7*M-1);
model.T3ALBRB = y(:,N+7*M:N+8*M-1);
model.TBGRB   = y(:,N+8*M:N+9*M-1);
model.TTRRB   = y(:,N+9*M:N+10*M-1);
model.ALBRB   = y(:,N+10*M:N+11*M-1);
model.T4RBT   = y(:,N+11*M:N+12*M-1);
model.T3RBT   = y(:,N+12*M:N+13*M-1);
model.XRB     = y(:,N+13*M:N+14*M-1);
model.XTTRRB  = y(:,N+14*M:N+15*M-1);
model.XTBGRB  = y(:,N+15*M:N+16*M-1);

model.fT4L    = y(:,N+16*M:N+17*M-1);
model.T4TBGL  = y(:,N+17*M:N+18*M-1);
model.T4TTRL  = y(:,N+18*M:N+19*M-1);
model.T4ALBL  = y(:,N+19*M:N+20*M-1);
model.fT3L    = y(:,N+20*M:N+21*M-1);
model.T3TBGL  = y(:,N+21*M:N+22*M-1);
model.T3TTRL  = y(:,N+22*M:N+23*M-1);
model.T3ALBL  = y(:,N+23*M:N+24*M-1);
model.TBGL    = y(:,N+24*M:N+25*M-1);
model.TTRL    = y(:,N+25*M:N+26*M-1);
model.ALBL    = y(:,N+26*M:N+27*M-1);
model.T4LT    = y(:,N+27*M:N+28*M-1);
model.T3LT    = y(:,N+28*M:N+29*M-1);
model.XL      = y(:,N+29*M:N+30*M-1);
model.XTTRL   = y(:,N+30*M:N+31*M-1);
model.XTBGL   = y(:,N+31*M:N+32*M-1);


%-----------Obtaining steady state values-----------%
% Body Blood Compartment
freeT4B = model.fT4B(end);
T4TBGB = model.T4TBGB(end);
T4TTRB = model.T4TTRB(end);
T4ALBB = model.T4ALBB(end);

freeT3B = model.fT3B(end);
T3TBGB = model.T3TBGB(end);
T3TTRB = model.T3TTRB(end);
T3ALBB = model.T3ALBB(end);

TBGB = model.TBGB(end);
TTRB = model.TTRB(end);
ALBB = model.ALBB(end);

% Thyroid blood Compartment
freeT4T = model.fT4T(end);
T4TBGT = model.T4TBGT(end);
T4TTRT = model.T4TTRT(end);
T4ALBT = model.T4ALBT(end);

freeT3T = model.fT3T(end);
T3TBGT = model.T3TBGT(end);
T3TTRT = model.T3TTRT(end);
T3ALBT = model.T3ALBT(end);

TBGT = model.TBGT(end);
TTRT = model.TTRT(end);
ALBT = model.ALBT(end);

% RB venous blood (last segment of RB blood) 
freeT4RB = model.fT4RB(end);
T4TBGRB = model.T4TBGRB(end);
T4TTRRB = model.T4TTRRB(end);
T4ALBRB = model.T4ALBRB(end);

freeT3RB = model.fT3RB(end);
T3TBGRB  = model.T3TBGRB(end);
T3TTRRB  =  model.T3TTRRB(end);
T3ALBRB  = model.T3ALBRB(end);

TBGRB = model.TBGRB(end);
TTRRB = model.TTRRB(end);
ALBRB = model.ALBRB(end);

% RB tissue (last segment)
T4RBT = model.T4RBT(end);
T3RBT = model.T3RBT(end);

% Liver venous blood (last segment of Liver blood) 
freeT4L = model.fT4L(end);
T4TBGL = model.T4TBGL(end);
T4TTRL = model.T4TTRL(end);
T4ALBL = model.T4ALBL(end);

freeT3L = model.fT3L(end);
T3TBGL = model.T3TBGL(end);
T3TTRL = model.T3TTRL(end);
T3ALBL = model.T3ALBL(end);

TBGL = model.TBGL(end);
TTRL = model.TTRL(end);
ALBL = model.ALBL(end);

%Liver tissue (last segment)
T4LT = model.T4LT(end);
T3LT = model.T3LT(end);



% Calculated metrics for T4 and T3 in Body Blood
boundT4B = T4TBGB + T4TTRB + T4ALBB;
totalT4B = boundT4B + freeT4B;
freeT4Bpercentage = freeT4B/(totalT4B)*100;
T4TBGBpercentage = T4TBGB/totalT4B*100;
T4TTRBpercentage = T4TTRB/totalT4B*100;
T4ALBBpercentage = T4ALBB/totalT4B*100;

TBGT4Bsatpercentage = T4TBGB/(TBGB+T4TBGB+T3TBGB)*100;
TTRT4Bsatpercentage = T4TTRB/(TTRB+T4TTRB+T3TTRB)*100;
ALBT4Bsatpercentage = T4ALBB/(ALBB+T4ALBB+T3ALBB)*100;

boundT3B = T3TBGB + T3TTRB + T3ALBB;
totalT3B = boundT3B + freeT3B;
freeT3Bpercentage = freeT3B/(totalT3B)*100;
T3TBGBpercentage = T3TBGB/totalT3B*100;
T3TTRBpercentage = T3TTRB/totalT3B*100;
T3ALBBpercentage = T3ALBB/totalT3B*100;

TBGT3Bsatpercentage = T3TBGB/(TBGB+T4TBGB+T3TBGB)*100;
TTRT3Bsatpercentage = T3TTRB/(TTRB+T4TTRB+T3TTRB)*100;
ALBT3Bsatpercentage = T3ALBB/(ALBB+T4ALBB+T3ALBB)*100;

TBGBsatpercentage = TBGT4Bsatpercentage + TBGT3Bsatpercentage;
TTRBsatpercentage = TTRT4Bsatpercentage + TTRT3Bsatpercentage;
ALBBsatpercentage = ALBT4Bsatpercentage + ALBT3Bsatpercentage;

% Calculated metrics for T4 and T3 in Thyroid blood
boundT4T = T4TBGT + T4TTRT + T4ALBT;
totalT4T = boundT4T + freeT4T;
freeT4Tpercentage = freeT4T/(totalT4T)*100;
T4TBGTpercentage = T4TBGT/totalT4T*100;
T4TTRTpercentage = T4TTRT/totalT4T*100;
T4ALBTpercentage = T4ALBT/totalT4T*100;

TBGT4Tsatpercentage = T4TBGT/(TBGT+T4TBGT+T3TBGT)*100;
TTRT4Tsatpercentage = T4TTRT/(TTRT+T4TTRT+T3TTRT)*100;
ALBT4Tsatpercentage = T4ALBT/(ALBT+T4ALBT+T3ALBT)*100;

boundT3T = T3TBGT + T3TTRT + T3ALBT;
totalT3T = boundT3T + freeT3T;
freeT3Tpercentage = freeT3T/(totalT3T)*100;
T3TBGTpercentage = T3TBGT/totalT3T*100;
T3TTRTpercentage = T3TTRT/totalT3T*100;
T3ALBTpercentage = T3ALBT/totalT3T*100;

TBGT3Tsatpercentage = T3TBGT/(TBGT+T4TBGT+T3TBGT)*100;
TTRT3Tsatpercentage = T3TTRT/(TTRT+T4TTRT+T3TTRT)*100;
ALBT3Tsatpercentage = T3ALBT/(ALBT+T4ALBT+T3ALBT)*100;

TBGTsatpercentage = TBGT4Tsatpercentage + TBGT3Tsatpercentage;
TTRTsatpercentage = TTRT4Tsatpercentage + TTRT3Tsatpercentage;
ALBTsatpercentage = ALBT4Tsatpercentage + ALBT3Tsatpercentage;

% Calculated metrics for T4 and T3 in RB venous blood (last segment of RB blood) 
boundT4RB = T4TBGRB + T4TTRRB + T4ALBRB;
totalT4RB = boundT4RB + freeT4RB;
freeT4RBpercentage = freeT4RB/(totalT4RB)*100;
T4TBGRBpercentage = T4TBGRB/totalT4RB*100;
T4TTRRBpercentage = T4TTRRB/totalT4RB*100;
T4ALBRBpercentage = T4ALBRB/totalT4RB*100;

TBGT4RBsatpercentage = T4TBGRB/(TBGRB+T4TBGRB+T3TBGRB)*100;
TTRT4RBsatpercentage = T4TTRRB/(TTRRB+T4TTRRB+T3TTRRB)*100;
ALBT4RBsatpercentage = T4ALBRB/(ALBRB+T4ALBRB+T3ALBRB)*100;

boundT3RB = T3TBGRB + T3TTRRB + T3ALBRB;
totalT3RB = boundT3RB + freeT3RB;
freeT3RBpercentage = freeT3RB/(totalT3RB)*100;
T3TBGRBpercentage = T3TBGRB/totalT3RB*100;
T3TTRRBpercentage = T3TTRRB/totalT3RB*100;
T3ALBRBpercentage = T3ALBRB/totalT3RB*100;

TBGT3RBsatpercentage = T3TBGRB/(TBGRB+T4TBGRB+T3TBGRB)*100;
TTRT3RBsatpercentage = T3TTRRB/(TTRRB+T4TTRRB+T3TTRRB)*100;
ALBT3RBsatpercentage = T3ALBRB/(ALBRB+T4ALBRB+T3ALBRB)*100;

TBGRBsatpercentage = TBGT4RBsatpercentage + TBGT3RBsatpercentage;
TTRRBsatpercentage = TTRT4RBsatpercentage + TTRT3RBsatpercentage;
ALBRBsatpercentage = ALBT4RBsatpercentage + ALBT3RBsatpercentage;

% Calculated metrics for T4 and T3 in Liver venous blood (last segment of Liver blood) 
boundT4L = T4TBGL + T4TTRL + T4ALBL;
totalT4L = boundT4L + freeT4L;
freeT4Lpercentage = freeT4L/(totalT4L)*100;
T4TBGLpercentage = T4TBGL/totalT4L*100;
T4TTRLpercentage = T4TTRL/totalT4L*100;
T4ALBLpercentage = T4ALBL/totalT4L*100;

TBGT4Lsatpercentage = T4TBGL/(TBGL+T4TBGL+T3TBGL)*100;
TTRT4Lsatpercentage = T4TTRL/(TTRL+T4TTRL+T3TTRL)*100;
ALBT4Lsatpercentage = T4ALBL/(ALBL+T4ALBL+T3ALBL)*100;

boundT3L = T3TBGL + T3TTRL + T3ALBL;
totalT3L = boundT3L + freeT3L;
freeT3Lpercentage = freeT3L/(totalT3L)*100;
T3TBGLpercentage = T3TBGL/totalT3L*100;
T3TTRLpercentage = T3TTRL/totalT3L*100;
T3ALBLpercentage = T3ALBL/totalT3L*100;

TBGT3Lsatpercentage = T3TBGL/(TBGL+T4TBGL+T3TBGL)*100;
TTRT3Lsatpercentage = T3TTRL/(TTRL+T4TTRL+T3TTRL)*100;
ALBT3Lsatpercentage = T3ALBL/(ALBL+T4ALBL+T3ALBL)*100;

TBGLsatpercentage = TBGT4Lsatpercentage + TBGT3Lsatpercentage;
TTRLsatpercentage = TTRT4Lsatpercentage + TTRT3Lsatpercentage;
ALBLsatpercentage = ALBT4Lsatpercentage + ALBT3Lsatpercentage;


% T4 and T3 differnetial between arterial and venous blood concentrations
CACVfT4T = freeT4T - freeT4B;
CACVfT4L = freeT4L - freeT4B;
CACVfT4RB = freeT4RB - freeT4B;

CACVfT4T_percent = (freeT4T - freeT4B)/freeT4B * 100;
CACVfT4L_percent = (freeT4L - freeT4B)/freeT4B * 100;
CACVfT4RB_percent = (freeT4RB - freeT4B)/freeT4B * 100;

CACVfT3T = freeT3T - freeT3B;
CACVfT3L = freeT3L - freeT3B;
CACVfT3RB = freeT3RB - freeT3B;

CACVfT3T_percent = (freeT3T - freeT3B)/freeT3B * 100;
CACVfT3L_percent = (freeT3L - freeT3B)/freeT3B * 100;
CACVfT3RB_percent = (freeT3RB - freeT3B)/freeT3B * 100;

CACVT4TBGT = T4TBGT - T4TBGB;
CACVT4TBGL = T4TBGL - T4TBGB;
CACVT4TBGRB = T4TBGRB - T4TBGB;

CACVT4TBGT_percent = (T4TBGT - T4TBGB)/T4TBGB * 100;
CACVT4TBGL_percent = (T4TBGL - T4TBGB)/T4TBGB * 100;
CACVT4TBGRB_percent = (T4TBGRB - T4TBGB)/T4TBGB * 100;

CACVT3TBGT = T3TBGT - T3TBGB;
CACVT3TBGL = T3TBGL - T3TBGB;
CACVT3TBGRB = T3TBGRB - T3TBGB;

CACVT3TBGT_percent = (T3TBGT - T3TBGB)/T3TBGB * 100;
CACVT3TBGL_percent = (T3TBGL - T3TBGB)/T3TBGB * 100;
CACVT3TBGRB_percent = (T3TBGRB - T3TBGB)/T3TBGB * 100;

CACVT4TTRT = T4TTRT - T4TTRB;
CACVT4TTRL = T4TTRL - T4TTRB;
CACVT4TTRRB = T4TTRRB - T4TTRB;

CACVT4TTRT_percent = (T4TTRT - T4TTRB)/T4TTRB * 100;
CACVT4TTRL_percent = (T4TTRL - T4TTRB)/T4TTRB * 100;
CACVT4TTRRB_percent = (T4TTRRB - T4TTRB)/T4TTRB * 100;

CACVT3TTRT = T3TTRT - T3TTRB;
CACVT3TTRL = T3TTRL - T3TTRB;
CACVT3TTRRB = T3TTRRB - T3TTRB;

CACVT3TTRT_percent = (T3TTRT - T3TTRB)/T3TTRB * 100;
CACVT3TTRL_percent = (T3TTRL - T3TTRB)/T3TTRB * 100;
CACVT3TTRRB_percent = (T3TTRRB - T3TTRB)/T3TTRB * 100;

CACVT4ALBT = T4ALBT - T4ALBB;
CACVT4ALBL = T4ALBL - T4ALBB;
CACVT4ALBRB = T4ALBRB - T4ALBB;

CACVT4ALBT_percent = (T4ALBT - T4ALBB)/T4ALBB * 100;
CACVT4ALBL_percent = (T4ALBL - T4ALBB)/T4ALBB * 100;
CACVT4ALBRB_percent = (T4ALBRB - T4ALBB)/T4ALBB * 100;

CACVT3ALBT = T3ALBT - T3ALBB;
CACVT3ALBL = T3ALBL - T3ALBB;
CACVT3ALBRB = T3ALBRB - T3ALBB;

CACVT3ALBT_percent = (T3ALBT - T3ALBB)/T3ALBB * 100;
CACVT3ALBL_percent = (T3ALBL - T3ALBB)/T3ALBB * 100;
CACVT3ALBRB_percent = (T3ALBRB - T3ALBB)/T3ALBB * 100;


% RB blood - all segments
freeT4RB = model.fT4RB(end, :);
T4TBGRB = model.T4TBGRB(end, :);
T4TTRRB = model.T4TTRRB(end, :);
T4ALBRB = model.T4ALBRB(end, :);

freeT3RB = model.fT3RB(end, :);
T3TBGRB  = model.T3TBGRB(end, :);
T3TTRRB  =  model.T3TTRRB(end, :);
T3ALBRB  = model.T3ALBRB(end, :);

TBGRB = model.TBGRB(end, :);
TTRRB = model.TTRRB(end, :);
ALBRB = model.ALBRB(end, :);

% RB tissue - all segments
T4RBT = model.T4RBT(end, :);
T3RBT = model.T3RBT(end, :);

% Liver blood - all segments
freeT4L = model.fT4L(end, :);
T4TBGL = model.T4TBGL(end, :);
T4TTRL = model.T4TTRL(end, :);
T4ALBL = model.T4ALBL(end, :);

freeT3L = model.fT3L(end, :);
T3TBGL = model.T3TBGL(end, :);
T3TTRL = model.T3TTRL(end, :);
T3ALBL = model.T3ALBL(end, :);

TBGL = model.TBGL(end, :);
TTRL = model.TTRL(end, :);
ALBL = model.ALBL(end, :);

% Liver tissue - all segments
T4LT = model.T4LT(end, :);
T3LT = model.T3LT(end, :);


% Binding Rates in Body Blood 
body_blood_T4TBG_assoc_rate = param.k1*freeT4B*TBGB;
body_blood_T4TBG_disassoc_rate = param.k2*T4TBGB;
body_blood_difference_T4TBG_rate = body_blood_T4TBG_disassoc_rate - body_blood_T4TBG_assoc_rate; 
body_blood_percent_difference_T4TBG_rate = body_blood_difference_T4TBG_rate/body_blood_T4TBG_assoc_rate*100;

body_blood_T4TTR_assoc_rate = param.k3*freeT4B*TTRB;
body_blood_T4TTR_disassoc_rate = param.k4*T4TTRB;
body_blood_difference_T4TTR_rate = body_blood_T4TTR_disassoc_rate - body_blood_T4TTR_assoc_rate; 
body_blood_percent_difference_T4TTR_rate = body_blood_difference_T4TTR_rate/body_blood_T4TTR_assoc_rate*100;

body_blood_T4ALB_assoc_rate = param.k5*freeT4B*ALBB;
body_blood_T4ALB_disassoc_rate = param.k6*T4ALBB;
body_blood_difference_T4ALB_rate = body_blood_T4ALB_disassoc_rate - body_blood_T4ALB_assoc_rate; 
body_blood_percent_difference_T4ALB_rate = body_blood_difference_T4ALB_rate/body_blood_T4ALB_assoc_rate*100;

body_blood_T3TBG_assoc_rate = param.k7*freeT3B*TBGB;
body_blood_T3TBG_disassoc_rate = param.k8*T3TBGB;
body_blood_difference_T3TBG_rate = body_blood_T3TBG_disassoc_rate - body_blood_T3TBG_assoc_rate; 
body_blood_percent_difference_T3TBG_rate = body_blood_difference_T3TBG_rate/body_blood_T3TBG_assoc_rate*100;

body_blood_T3TTR_assoc_rate = param.k9*freeT3B*TTRB;
body_blood_T3TTR_disassoc_rate = param.k10*T3TTRB;
body_blood_difference_T3TTR_rate = body_blood_T3TTR_disassoc_rate - body_blood_T3TTR_assoc_rate; 
body_blood_percent_difference_T3TTR_rate = body_blood_difference_T3TTR_rate/body_blood_T3TTR_assoc_rate*100;

body_blood_T3ALB_assoc_rate = param.k11*freeT3B*ALBB;
body_blood_T3ALB_disassoc_rate = param.k12*T3ALBB;
body_blood_difference_T3ALB_rate = body_blood_T3ALB_disassoc_rate - body_blood_T3ALB_assoc_rate; 
body_blood_percent_difference_T3ALB_rate = body_blood_difference_T3ALB_rate/body_blood_T3ALB_assoc_rate*100;

% Binding Rates in Thyroid blood 
thyroid_T4TBG_assoc_rate = param.k1*freeT4T*TBGT;
thyroid_T4TBG_disassoc_rate = param.k2*T4TBGT;
thyroid_difference_T4TBG_rate = thyroid_T4TBG_disassoc_rate - thyroid_T4TBG_assoc_rate; 
thyroid_percent_difference_T4TBG_rate = thyroid_difference_T4TBG_rate/thyroid_T4TBG_assoc_rate*100;

thyroid_T4TTR_assoc_rate = param.k3*freeT4T*TTRT;
thyroid_T4TTR_disassoc_rate = param.k4*T4TTRT;
thyroid_difference_T4TTR_rate = thyroid_T4TTR_disassoc_rate - thyroid_T4TTR_assoc_rate; 
thyroid_percent_difference_T4TTR_rate = thyroid_difference_T4TTR_rate/thyroid_T4TTR_assoc_rate*100;

thyroid_T4ALB_assoc_rate = param.k5*freeT4T*ALBT;
thyroid_T4ALB_disassoc_rate = param.k6*T4ALBT;
thyroid_difference_T4ALB_rate = thyroid_T4ALB_disassoc_rate - thyroid_T4ALB_assoc_rate; 
thyroid_percent_difference_T4ALB_rate = thyroid_difference_T4ALB_rate/thyroid_T4ALB_assoc_rate*100;

thyroid_T3TBG_assoc_rate = param.k7*freeT3T*TBGT;
thyroid_T3TBG_disassoc_rate = param.k8*T3TBGT;
thyroid_difference_T3TBG_rate = thyroid_T3TBG_disassoc_rate - thyroid_T3TBG_assoc_rate; 
thyroid_percent_difference_T3TBG_rate = thyroid_difference_T3TBG_rate/thyroid_T3TBG_assoc_rate*100;

thyroid_T3TTR_assoc_rate = param.k9*freeT3T*TTRT;
thyroid_T3TTR_disassoc_rate = param.k10*T3TTRT;
thyroid_difference_T3TTR_rate = thyroid_T3TTR_disassoc_rate - thyroid_T3TTR_assoc_rate; 
thyroid_percent_difference_T3TTR_rate = thyroid_difference_T3TTR_rate/thyroid_T3TTR_assoc_rate*100;

thyroid_T3ALB_assoc_rate = param.k11*freeT3T*ALBT;
thyroid_T3ALB_disassoc_rate = param.k12*T3ALBT;
thyroid_difference_T3ALB_rate = thyroid_T3ALB_disassoc_rate - thyroid_T3ALB_assoc_rate; 
thyroid_percent_difference_T3ALB_rate = thyroid_difference_T3ALB_rate/thyroid_T3ALB_assoc_rate*100;

% Binding Rates in Liver blood 
liver_T4TBG_assoc_rate = param.k1.*freeT4L.*TBGL;
liver_T4TBG_disassoc_rate = param.k2.*T4TBGL;
liver_difference_T4TBG_rate = liver_T4TBG_disassoc_rate - liver_T4TBG_assoc_rate; 
liver_percent_difference_T4TBG_rate = liver_difference_T4TBG_rate./liver_T4TBG_assoc_rate*100;

liver_T4TTR_assoc_rate = param.k3.*freeT4L.*TTRL;
liver_T4TTR_disassoc_rate = param.k4.*T4TTRL;
liver_difference_T4TTR_rate = liver_T4TTR_disassoc_rate - liver_T4TTR_assoc_rate; 
liver_percent_difference_T4TTR_rate = liver_difference_T4TTR_rate./liver_T4TTR_assoc_rate*100;

liver_T4ALB_assoc_rate = param.k5.*freeT4L.*ALBL;
liver_T4ALB_disassoc_rate = param.k6.*T4ALBL;
liver_difference_T4ALB_rate = liver_T4ALB_disassoc_rate - liver_T4ALB_assoc_rate; 
liver_percent_difference_T4ALB_rate = liver_difference_T4ALB_rate./liver_T4ALB_assoc_rate*100;

liver_T3TBG_assoc_rate = param.k7.*freeT3L.*TBGL;
liver_T3TBG_disassoc_rate = param.k8.*T3TBGL;
liver_difference_T3TBG_rate = liver_T3TBG_disassoc_rate - liver_T3TBG_assoc_rate; 
liver_percent_difference_T3TBG_rate = liver_difference_T3TBG_rate./liver_T3TBG_assoc_rate*100;

liver_T3TTR_assoc_rate = param.k9.*freeT3L.*TTRL;
liver_T3TTR_disassoc_rate = param.k10.*T3TTRL;
liver_difference_T3TTR_rate = liver_T3TTR_disassoc_rate - liver_T3TTR_assoc_rate; 
liver_percent_difference_T3TTR_rate = liver_difference_T3TTR_rate./liver_T3TTR_assoc_rate*100;

liver_T3ALB_assoc_rate = param.k11.*freeT3L.*ALBL;
liver_T3ALB_disassoc_rate = param.k12.*T3ALBL;
liver_difference_T3ALB_rate = liver_T3ALB_disassoc_rate - liver_T3ALB_assoc_rate; 
liver_percent_difference_T3ALB_rate = liver_difference_T3ALB_rate./liver_T3ALB_assoc_rate*100;

% Binding Rates in RB blood 
RB_T4TBG_assoc_rate = param.k1.*freeT4RB.*TBGRB;
RB_T4TBG_disassoc_rate = param.k2.*T4TBGRB;
RB_difference_T4TBG_rate = RB_T4TBG_disassoc_rate - RB_T4TBG_assoc_rate; 
RB_percent_difference_T4TBG_rate = RB_difference_T4TBG_rate./RB_T4TBG_assoc_rate*100;

RB_T4TTR_assoc_rate = param.k3.*freeT4RB.*TTRRB;
RB_T4TTR_disassoc_rate = param.k4.*T4TTRRB;
RB_difference_T4TTR_rate = RB_T4TTR_disassoc_rate - RB_T4TTR_assoc_rate; 
RB_percent_difference_T4TTR_rate = RB_difference_T4TTR_rate./RB_T4TTR_assoc_rate*100;

RB_T4ALB_assoc_rate = param.k5.*freeT4RB.*ALBRB;
RB_T4ALB_disassoc_rate = param.k6.*T4ALBRB;
RB_difference_T4ALB_rate = RB_T4ALB_disassoc_rate - RB_T4ALB_assoc_rate; 
RB_percent_difference_T4ALB_rate = RB_difference_T4ALB_rate./RB_T4ALB_assoc_rate*100;

RB_T3TBG_assoc_rate = param.k7.*freeT3RB.*TBGRB;
RB_T3TBG_disassoc_rate = param.k8.*T3TBGRB;
RB_difference_T3TBG_rate = RB_T3TBG_disassoc_rate - RB_T3TBG_assoc_rate; 
RB_percent_difference_T3TBG_rate = RB_difference_T3TBG_rate./RB_T3TBG_assoc_rate*100;

RB_T3TTR_assoc_rate = param.k9.*freeT3RB.*TTRRB;
RB_T3TTR_disassoc_rate = param.k10.*T3TTRRB;
RB_difference_T3TTR_rate = RB_T3TTR_disassoc_rate - RB_T3TTR_assoc_rate; 
RB_percent_difference_T3TTR_rate = RB_difference_T3TTR_rate./RB_T3TTR_assoc_rate*100;

RB_T3ALB_assoc_rate = param.k11.*freeT3RB.*ALBRB;
RB_T3ALB_disassoc_rate = param.k12.*T3ALBRB;
RB_difference_T3ALB_rate = RB_T3ALB_disassoc_rate - RB_T3ALB_assoc_rate; 
RB_percent_difference_T3ALB_rate = RB_difference_T3ALB_rate./RB_T3ALB_assoc_rate*100;


%% ------------------------- Plotting Results (Fig. 3 and dashed lines in Fig. 7) ------- %%

%Free T4 in Liver and RB blood (Fig. 3A)
figure(30)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("fT4")

%Free T4 in Liver blood (Fig. 7A dashed line)
figure(71)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present")
title("Liver blood fT4")

%Free T4 in RB blood (Fig. 7B dashed line)
figure(72)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present")
title("RB blood fT4")

%T4TBG in Liver and RB blood (Fig. 3B)
figure(31)
hold on
plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("T4TBG")

%T4TTR in Liver and RB blood (Fig. 3C)
figure(32)
hold on
plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("T4TTR")

%T4ALB in Liver and RB blood (Fig. 3D)
figure(33)
hold on
plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("T4ALB")

%Total T4 in Liver and RB blood (Fig. 3E)
figure(34)
hold on
plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4L(end,:)+model.T4TBGL(end,:)+model.T4TTRL(end,:)+model.T4ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4RB(end,:)+model.T4TBGRB(end,:)+model.T4TTRRB(end,:)+model.T4ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("Total T4")

%Free T3 in Liver and RB blood (Fig. 3F)
figure(35)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("fT3")

%Free T3 in Liver blood (Fig. 7E dashed line)
figure(75)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present")
title("Liver blood fT3")

%Free T3 in RB blood (Fig. 7F dashed line)
figure(76)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present")
title("RB blood fT3")

%T3TBG in Liver and RB blood (Fig. 3G)
figure(36)
hold on
plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("T3TBG")

%T3TTR in Liver and RB blood (Fig. 3H)
figure(37)
hold on
plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("T3TTR")

%T3ALB in Liver and RB blood (Fig. 3I)
figure(38)
hold on
plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("T3ALB")

%Total T3 in Liver and RB blood (Fig. 3J)
figure(39)
hold on
plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3L(end,:)+model.T3TBGL(end,:)+model.T3TTRL(end,:)+model.T3ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3RB(end,:)+model.T3TBGRB(end,:)+model.T3TTRRB(end,:)+model.T3ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("Total T3")

%T4 in Liver tissue  (Fig. 7C dashed line)
figure(73)
hold on
plot([1:1:param.M],model.T4LT(end,:), 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present")
title("Liver tissue T4")

%T4 in RB tissue  (Fig. 7D dashed line)
figure(74)
hold on
plot([1:1:param.M],model.T4RBT(end,:), 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present")
title("RB tissue T4")

%T3 in Liver tissue (Fig. 7G dashed line)
figure(77)
hold on
plot([1:1:param.M],model.T3LT(end,:),'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present")
title("Liver tissue T3")

%T3 in RB tissue (Fig. 7H dashed line)
figure(78)
hold on
plot([1:1:param.M],model.T3RBT(end,:), 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present")
title("RB tissue T3")

%Free TBG in Liver and RB Blood (not used)
figure(305)
hold on
plot([0:1:param.M],[model.TBGB(end),model.TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1,'LineStyle', '--')
plot([0:1:param.M],[model.TBGB(end),model.TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1,'LineStyle', '--')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("Free TBG")

%Free TTR in Liver and RB Blood (not used)
figure(306)
hold on
plot([0:1:param.M],[model.TTRB(end),model.TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.TTRB(end),model.TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("Free TTR")

%Free ALB in Liver and RB Blood (not used)
figure(307)
hold on
plot([0:1:param.M],[model.ALBB(end),model.ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1)
plot([0:1:param.M],[model.ALBB(end),model.ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1)
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver blood", "RB blood")
title("Free ALB")

%% ------------------------- Plotting Results (Fig. 4) ---------------------------------- %%

%Absolute Drop of T4 in Liver blood (Fig. 4A)
figure(41)
hold on
plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBL(end,1) + model.T4TTRB(end)-model.T4TTRL(end,1) + model.T4TBGB(end)-model.T4TBGL(end,1) + model.fT4B(end)-model.fT4L(end,1), - model.T4TBGL(end,2:param.M) + model.T4TBGL(end,1:param.M-1) - model.T4ALBL(end,2:param.M) + model.T4ALBL(end,1:param.M-1) - model.T4TTRL(end,2:param.M) + model.T4TTRL(end,1:param.M-1) - model.fT4L(end,2:param.M) + model.fT4L(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBL(end,1), - model.T4ALBL(end,2:param.M) + model.T4ALBL(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T4TTRB(end)-model.T4TTRL(end,1), - model.T4TTRL(end,2:param.M) + model.T4TTRL(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T4TBGB(end)-model.T4TBGL(end,1), - model.T4TBGL(end,2:param.M) + model.T4TBGL(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
plot([1:1:param.M], [model.fT4B(end)-model.fT4L(end,1), - model.fT4L(end,2:param.M) + model.fT4L(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
xlabel('Tissue Segment')
ylabel('Plasma Concentration Drop (pM)')
legend("Total T4", "T4ALB", "T4TTR", "T4TBG", "fT4")
title("Liver blood T4")

%Absolute Drop of T4 in RB blood (Fig. 4B)
figure(42)
hold on
plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBRB(end,1) + model.T4TTRB(end)-model.T4TTRRB(end,1) + model.T4TBGB(end)-model.T4TBGRB(end,1) + model.fT4B(end)-model.fT4RB(end,1), - model.T4TBGRB(end,2:param.M) + model.T4TBGRB(end,1:param.M-1) - model.T4ALBRB(end,2:param.M) + model.T4ALBRB(end,1:param.M-1) - model.T4TTRRB(end,2:param.M) + model.T4TTRRB(end,1:param.M-1) - model.fT4RB(end,2:param.M) + model.fT4RB(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBRB(end,1), - model.T4ALBRB(end,2:param.M) + model.T4ALBRB(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T4TTRB(end)-model.T4TTRRB(end,1), - model.T4TTRRB(end,2:param.M) + model.T4TTRRB(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T4TBGB(end)-model.T4TBGRB(end,1), - model.T4TBGRB(end,2:param.M) + model.T4TBGRB(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
plot([1:1:param.M], [model.fT4B(end)-model.fT4RB(end,1), - model.fT4RB(end,2:param.M) + model.fT4RB(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
xlabel('Tissue Segment')
ylabel('Plasma Concentration Drop (pM)')
legend("Total T4", "T4ALB", "T4TTR", "T4TBG", "fT4")
title("RB blood T4")

%Absolute Drop of T3 in Liver blood (Fig. 4C)
figure(43)
hold on
plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBL(end,1) + model.T3TTRB(end)-model.T3TTRL(end,1) + model.T3TBGB(end)-model.T3TBGL(end,1) + model.fT3B(end)-model.fT3L(end,1), - model.T3TBGL(end,2:param.M) + model.T3TBGL(end,1:param.M-1) - model.T3ALBL(end,2:param.M) + model.T3ALBL(end,1:param.M-1) - model.T3TTRL(end,2:param.M) + model.T3TTRL(end,1:param.M-1) - model.fT3L(end,2:param.M) + model.fT3L(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBL(end,1), - model.T3ALBL(end,2:param.M) + model.T3ALBL(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T3TTRB(end)-model.T3TTRL(end,1), - model.T3TTRL(end,2:param.M) + model.T3TTRL(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T3TBGB(end)-model.T3TBGL(end,1), - model.T3TBGL(end,2:param.M) + model.T3TBGL(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
plot([1:1:param.M], [model.fT3B(end)-model.fT3L(end,1), - model.fT3L(end,2:param.M) + model.fT3L(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
xlabel('Tissue Segment')
ylabel('Plasma Concentration Drop (pM)')
legend("Total T3", "T3ALB", "T3TTR", "T3TBG", "fT3")
title("Liver blood T3")

%Absolute Drop of T3 in RB blood (Fig. 4D)
figure(44)
hold on
plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBRB(end,1) + model.T3TTRB(end)-model.T3TTRRB(end,1) + model.T3TBGB(end)-model.T3TBGRB(end,1) + model.fT3B(end)-model.fT3RB(end,1), - model.T3TBGRB(end,2:param.M) + model.T3TBGRB(end,1:param.M-1) - model.T3ALBRB(end,2:param.M) + model.T3ALBRB(end,1:param.M-1) - model.T3TTRRB(end,2:param.M) + model.T3TTRRB(end,1:param.M-1) - model.fT3RB(end,2:param.M) + model.fT3RB(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBRB(end,1), - model.T3ALBRB(end,2:param.M) + model.T3ALBRB(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T3TTRB(end)-model.T3TTRRB(end,1), - model.T3TTRRB(end,2:param.M) + model.T3TTRRB(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
plot([1:1:param.M], [model.T3TBGB(end)-model.T3TBGRB(end,1), - model.T3TBGRB(end,2:param.M) + model.T3TBGRB(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
plot([1:1:param.M], [model.fT3B(end)-model.fT3RB(end,1), - model.fT3RB(end,2:param.M) + model.fT3RB(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
xlabel('Tissue Segment')
ylabel('Plasma Concentration Drop (pM)')
legend("Total T3", "T3ALB", "T3TTR", "T3TBG", "fT3")
title("RB blood T3")

%% ------------------------- Plotting Results (Fig. 5) ---------------------------------- %%

%Liver blood T4TBG (Fig. 5A)
figure(51)
hold on
yyaxis left
plot([0:1:param.M],[body_blood_T4TBG_assoc_rate, liver_T4TBG_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T4TBG_disassoc_rate, liver_T4TBG_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T4TBG_rate, liver_percent_difference_T4TBG_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Distance from Equilibrium (%)')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT4+TBG<->T4TBG')


%Liver blood T4TTR (Fig. 5B)
figure(52)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T4TTR_assoc_rate, liver_T4TTR_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T4TTR_disassoc_rate, liver_T4TTR_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T4TTR_rate, liver_percent_difference_T4TTR_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Distance from Equilibrium (%)')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT4+TTR<->T4TTR')

%Liver blood T4ALB (Fig. 5C)
figure(53)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T4ALB_assoc_rate, liver_T4ALB_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T4ALB_disassoc_rate, liver_T4ALB_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T4ALB_rate liver_percent_difference_T4ALB_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Distance from Equilibrium (%)')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT4+ALB<->T4ALB')

%RB blood T4TBG (Fig. 5D)
figure(54)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T4TBG_assoc_rate, RB_T4TBG_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T4TBG_disassoc_rate, RB_T4TBG_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T4TBG_rate, RB_percent_difference_T4TBG_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Distance from Equilibrium (%)')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT4+TBG<->T4TBG')

%RB blood T4TTR (Fig. 5E)
figure(55)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T4TTR_assoc_rate, RB_T4TTR_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T4TTR_disassoc_rate, RB_T4TTR_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T4TTR_rate, RB_percent_difference_T4TTR_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Distance from Equilibrium (%)')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT4+TTR<->T4TTR')

%RB blood T4ALB (Fig. 5F)
figure(56)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T4ALB_assoc_rate, RB_T4ALB_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T4ALB_disassoc_rate, RB_T4ALB_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Association or Dissociation Rate (pM/S)')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T4ALB_rate, RB_percent_difference_T4ALB_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Distance from Equilibrium (%)')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT4+ALB<->T4ALB')

%% ------------------------- Plotting Results (Fig. 6) ---------------------------------- %%

%Liver blood T3TBG (Fig. 6A)
figure(61)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T3TBG_assoc_rate, liver_T3TBG_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T3TBG_disassoc_rate, liver_T3TBG_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T3TBG_rate, liver_percent_difference_T3TBG_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Percent Difference')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT3+TBG<->T3TBG')

%Liver blood T3TTR (Fig. 6B)
figure(62)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T3TTR_assoc_rate, liver_T3TTR_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T3TTR_disassoc_rate, liver_T3TTR_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T3TTR_rate, liver_percent_difference_T3TTR_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Percent Difference')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT3+TTR<->T3TTR')

%Liver blood T3ALB (Fig. 6C)
figure(63)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T3ALB_assoc_rate, liver_T3ALB_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T3ALB_disassoc_rate, liver_T3ALB_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T3ALB_rate, liver_percent_difference_T3ALB_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Percent Difference')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT3+ALB<->T3ALB')

%RB blood T3TBG (Fig. 6D)
figure(64)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T3TBG_assoc_rate, RB_T3TBG_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T3TBG_disassoc_rate, RB_T3TBG_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T3TBG_rate, RB_percent_difference_T3TBG_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Percent Difference')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT3+TBG<->T3TBG')

%RB blood T3TTR (Fig. 6E)
figure(65)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T3TTR_assoc_rate, RB_T3TTR_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T3TTR_disassoc_rate, RB_T3TTR_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T3TTR_rate, RB_percent_difference_T3TTR_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Percent Difference')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT3+TTR<->T3TTR')

%RB blood T3ALB (Fig. 6F)
figure(66)
hold on
yyaxis left
plot([0:1:param.M], [body_blood_T3ALB_assoc_rate, RB_T3ALB_assoc_rate], 'Color', [0,0.45,0.74], 'LineWidth', 1.5, 'LineStyle', '-')
plot([0:1:param.M], [body_blood_T3ALB_disassoc_rate, RB_T3ALB_disassoc_rate], 'Color', [0.47,0.67,0.19], 'LineWidth', 1.5, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Association or Dissociation Rate (pM/S)')
yyaxis right
plot([0:1:param.M], [body_blood_percent_difference_T3ALB_rate, RB_percent_difference_T3ALB_rate], 'Color', [0.85,0.33,0.1], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('Percent Difference')
legend('Association rate', 'Dissociation rate', 'Distance from equilibrium')
title('fT3+ALB<->T3ALB')

%% ------------------------- Plotting Results (Solid lines in Fig. 7) ------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Setting all THBPs to zero
param.TBGtot = 0;  
param.TTRtot = 0;   
param.ALBtot = 0;

% Update the initial THBPs concentrations accordingly
init.TBGB    = param.TBGtot/param.VB;
init.TTRB    = param.TTRtot/param.VB;
init.ALBB    = param.ALBtot/param.VB;

y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600];  %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode15s('TH_PBK_Spatial_ODE',tspan, y0, options, param); 

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.XB      = y(:,23);
model.XTTRB   = y(:,24);
model.XT      = y(:,25);
model.XTTRT   = y(:,26);
model.XTBGB   = y(:,27);
model.XTBGT   = y(:,28);

N = 29;  %Starting index of repeating variables
M = param.M;
model.fT4RB   = y(:,N:N+M-1);
model.T4TBGRB = y(:,N+M:N+2*M-1);
model.T4TTRRB = y(:,N+2*M:N+3*M-1);
model.T4ALBRB = y(:,N+3*M:N+4*M-1);
model.fT3RB   = y(:,N+4*M:N+5*M-1);
model.T3TBGRB = y(:,N+5*M:N+6*M-1);
model.T3TTRRB = y(:,N+6*M:N+7*M-1);
model.T3ALBRB = y(:,N+7*M:N+8*M-1);
model.TBGRB   = y(:,N+8*M:N+9*M-1);
model.TTRRB   = y(:,N+9*M:N+10*M-1);
model.ALBRB   = y(:,N+10*M:N+11*M-1);
model.T4RBT   = y(:,N+11*M:N+12*M-1);
model.T3RBT   = y(:,N+12*M:N+13*M-1);
model.XRB     = y(:,N+13*M:N+14*M-1);
model.XTTRRB  = y(:,N+14*M:N+15*M-1);
model.XTBGRB  = y(:,N+15*M:N+16*M-1);

model.fT4L    = y(:,N+16*M:N+17*M-1);
model.T4TBGL  = y(:,N+17*M:N+18*M-1);
model.T4TTRL  = y(:,N+18*M:N+19*M-1);
model.T4ALBL  = y(:,N+19*M:N+20*M-1);
model.fT3L    = y(:,N+20*M:N+21*M-1);
model.T3TBGL  = y(:,N+21*M:N+22*M-1);
model.T3TTRL  = y(:,N+22*M:N+23*M-1);
model.T3ALBL  = y(:,N+23*M:N+24*M-1);
model.TBGL    = y(:,N+24*M:N+25*M-1);
model.TTRL    = y(:,N+25*M:N+26*M-1);
model.ALBL    = y(:,N+26*M:N+27*M-1);
model.T4LT    = y(:,N+27*M:N+28*M-1);
model.T3LT    = y(:,N+28*M:N+29*M-1);
model.XL      = y(:,N+29*M:N+30*M-1);
model.XTTRL   = y(:,N+30*M:N+31*M-1);
model.XTBGL   = y(:,N+31*M:N+32*M-1);

%Free T4 in Liver blood (Fig. 7A solid line)
figure(71)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present","THBPs Absent")
title("Liver blood fT4")

%Free T4 in RB blood (Fig. 7B solid line)
figure(72)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present","THBPs Absent")
title("RB blood fT4")

%T4 in Liver tissue  (Fig. 7C solid line)
figure(73)
hold on
plot([1:1:param.M],model.T4LT(end,:), 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present","THBPs Absent")
title("Liver tissue T4")

%T4 in RB tissue  (Fig. 7D solid line)
figure(74)
hold on
plot([1:1:param.M],model.T4RBT(end,:), 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present","THBPs Absent")
title("RB tissue T4")

%Free T3 in Liver blood (Fig. 7E solid line)
figure(75)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present","THBPs Absent")
title("Liver blood fT3")

%Free T3 in RB blood (Fig. 7F solid line)
figure(76)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present","THBPs Absent")
title("RB blood fT3")

%T3 in Liver tissue (Fig. 7G solid line)
figure(77)
hold on
plot([1:1:param.M],model.T3LT(end,:),'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present","THBPs Absent")
title("Liver tissue T3")

%T3 in RB tissue (Fig. 7H solid line)
figure(78)
hold on
plot([1:1:param.M],model.T3RBT(end,:), 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("THBPs Present","THBPs Absent")
title("RB tissue T3")

%% ------------------------- Plotting Results (Fig. S2) --------------------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Setting TTR and ALB to zero 
param.TTRtot = 0;   
param.ALBtot = 0;

% Update the initial THBPs concentrations accordingly
init.TTRB    = param.TTRtot/param.VB;
init.ALBB    = param.ALBtot/param.VB;

y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600];  %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode15s('TH_PBK_Spatial_ODE',tspan, y0, options, param); 

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.XB      = y(:,23);
model.XTTRB   = y(:,24);
model.XT      = y(:,25);
model.XTTRT   = y(:,26);
model.XTBGB   = y(:,27);
model.XTBGT   = y(:,28);

N = 29;  %Starting index of repeating variables
M = param.M;
model.fT4RB   = y(:,N:N+M-1);
model.T4TBGRB = y(:,N+M:N+2*M-1);
model.T4TTRRB = y(:,N+2*M:N+3*M-1);
model.T4ALBRB = y(:,N+3*M:N+4*M-1);
model.fT3RB   = y(:,N+4*M:N+5*M-1);
model.T3TBGRB = y(:,N+5*M:N+6*M-1);
model.T3TTRRB = y(:,N+6*M:N+7*M-1);
model.T3ALBRB = y(:,N+7*M:N+8*M-1);
model.TBGRB   = y(:,N+8*M:N+9*M-1);
model.TTRRB   = y(:,N+9*M:N+10*M-1);
model.ALBRB   = y(:,N+10*M:N+11*M-1);
model.T4RBT   = y(:,N+11*M:N+12*M-1);
model.T3RBT   = y(:,N+12*M:N+13*M-1);
model.XRB     = y(:,N+13*M:N+14*M-1);
model.XTTRRB  = y(:,N+14*M:N+15*M-1);
model.XTBGRB  = y(:,N+15*M:N+16*M-1);

model.fT4L    = y(:,N+16*M:N+17*M-1);
model.T4TBGL  = y(:,N+17*M:N+18*M-1);
model.T4TTRL  = y(:,N+18*M:N+19*M-1);
model.T4ALBL  = y(:,N+19*M:N+20*M-1);
model.fT3L    = y(:,N+20*M:N+21*M-1);
model.T3TBGL  = y(:,N+21*M:N+22*M-1);
model.T3TTRL  = y(:,N+22*M:N+23*M-1);
model.T3ALBL  = y(:,N+23*M:N+24*M-1);
model.TBGL    = y(:,N+24*M:N+25*M-1);
model.TTRL    = y(:,N+25*M:N+26*M-1);
model.ALBL    = y(:,N+26*M:N+27*M-1);
model.T4LT    = y(:,N+27*M:N+28*M-1);
model.T3LT    = y(:,N+28*M:N+29*M-1);
model.XL      = y(:,N+29*M:N+30*M-1);
model.XTTRL   = y(:,N+30*M:N+31*M-1);
model.XTBGL   = y(:,N+31*M:N+32*M-1);


%Free T4 in Liver and RB blood (Fig. S2A)
figure(201)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT4 in tissue blood with TBG only")

%T4TBG in Liver and RB blood (Fig. S2B)
figure(202)
hold on
plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4TBG in tissue blood with TBG only")

%Free T3 in Liver and RB blood (Fig. S2C)
figure(203)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT3 in tissue blood with TBG only")

%T3TBG in Liver and RB blood (Fig. S2D)
figure(204)
hold on
plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3TBG in tissue blood with TBG only")


%% ------------------------- Plotting Results (Fig. S3) --------------------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Setting TBG and ALB to zero 
param.TBGtot = 0;   
param.ALBtot = 0;

% Update the initial THBPs concentrations accordingly
init.TBGB    = param.TBGtot/param.VB;
init.ALBB    = param.ALBtot/param.VB;

y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600];  %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode15s('TH_PBK_Spatial_ODE',tspan, y0, options, param); 

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.XB      = y(:,23);
model.XTTRB   = y(:,24);
model.XT      = y(:,25);
model.XTTRT   = y(:,26);
model.XTBGB   = y(:,27);
model.XTBGT   = y(:,28);

N = 29;  %Starting index of repeating variables
M = param.M;
model.fT4RB   = y(:,N:N+M-1);
model.T4TBGRB = y(:,N+M:N+2*M-1);
model.T4TTRRB = y(:,N+2*M:N+3*M-1);
model.T4ALBRB = y(:,N+3*M:N+4*M-1);
model.fT3RB   = y(:,N+4*M:N+5*M-1);
model.T3TBGRB = y(:,N+5*M:N+6*M-1);
model.T3TTRRB = y(:,N+6*M:N+7*M-1);
model.T3ALBRB = y(:,N+7*M:N+8*M-1);
model.TBGRB   = y(:,N+8*M:N+9*M-1);
model.TTRRB   = y(:,N+9*M:N+10*M-1);
model.ALBRB   = y(:,N+10*M:N+11*M-1);
model.T4RBT   = y(:,N+11*M:N+12*M-1);
model.T3RBT   = y(:,N+12*M:N+13*M-1);
model.XRB     = y(:,N+13*M:N+14*M-1);
model.XTTRRB  = y(:,N+14*M:N+15*M-1);
model.XTBGRB  = y(:,N+15*M:N+16*M-1);

model.fT4L    = y(:,N+16*M:N+17*M-1);
model.T4TBGL  = y(:,N+17*M:N+18*M-1);
model.T4TTRL  = y(:,N+18*M:N+19*M-1);
model.T4ALBL  = y(:,N+19*M:N+20*M-1);
model.fT3L    = y(:,N+20*M:N+21*M-1);
model.T3TBGL  = y(:,N+21*M:N+22*M-1);
model.T3TTRL  = y(:,N+22*M:N+23*M-1);
model.T3ALBL  = y(:,N+23*M:N+24*M-1);
model.TBGL    = y(:,N+24*M:N+25*M-1);
model.TTRL    = y(:,N+25*M:N+26*M-1);
model.ALBL    = y(:,N+26*M:N+27*M-1);
model.T4LT    = y(:,N+27*M:N+28*M-1);
model.T3LT    = y(:,N+28*M:N+29*M-1);
model.XL      = y(:,N+29*M:N+30*M-1);
model.XTTRL   = y(:,N+30*M:N+31*M-1);
model.XTBGL   = y(:,N+31*M:N+32*M-1);


%Free T4 in Liver and RB blood (Fig. S3A)
figure(301)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT4 in tissue blood with TTR only")

%T4TTR in Liver and RB blood (Fig. S3B)
figure(302)
hold on
plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4TTR in tissue blood with TTR only")

%Free T3 in Liver and RB blood (Fig. S3C)
figure(303)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT3 in tissue blood with TTR only")

%T3TTR in Liver and RB blood (Fig. S3D)
figure(304)
hold on
plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3TTR in tissue blood with TTR only")


%% ------------------------- Plotting Results (Fig. S4) --------------------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Setting TBG and TTR to zero 
param.TBGtot = 0;   
param.TTRtot = 0;

% Update the initial THBPs concentrations accordingly
init.TBGB    = param.TBGtot/param.VB;
init.TTRB    = param.TTRtot/param.VB;

y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600];  %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode15s('TH_PBK_Spatial_ODE',tspan, y0, options, param); 

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.XB      = y(:,23);
model.XTTRB   = y(:,24);
model.XT      = y(:,25);
model.XTTRT   = y(:,26);
model.XTBGB   = y(:,27);
model.XTBGT   = y(:,28);

N = 29;  %Starting index of repeating variables
M = param.M;
model.fT4RB   = y(:,N:N+M-1);
model.T4TBGRB = y(:,N+M:N+2*M-1);
model.T4TTRRB = y(:,N+2*M:N+3*M-1);
model.T4ALBRB = y(:,N+3*M:N+4*M-1);
model.fT3RB   = y(:,N+4*M:N+5*M-1);
model.T3TBGRB = y(:,N+5*M:N+6*M-1);
model.T3TTRRB = y(:,N+6*M:N+7*M-1);
model.T3ALBRB = y(:,N+7*M:N+8*M-1);
model.TBGRB   = y(:,N+8*M:N+9*M-1);
model.TTRRB   = y(:,N+9*M:N+10*M-1);
model.ALBRB   = y(:,N+10*M:N+11*M-1);
model.T4RBT   = y(:,N+11*M:N+12*M-1);
model.T3RBT   = y(:,N+12*M:N+13*M-1);
model.XRB     = y(:,N+13*M:N+14*M-1);
model.XTTRRB  = y(:,N+14*M:N+15*M-1);
model.XTBGRB  = y(:,N+15*M:N+16*M-1);

model.fT4L    = y(:,N+16*M:N+17*M-1);
model.T4TBGL  = y(:,N+17*M:N+18*M-1);
model.T4TTRL  = y(:,N+18*M:N+19*M-1);
model.T4ALBL  = y(:,N+19*M:N+20*M-1);
model.fT3L    = y(:,N+20*M:N+21*M-1);
model.T3TBGL  = y(:,N+21*M:N+22*M-1);
model.T3TTRL  = y(:,N+22*M:N+23*M-1);
model.T3ALBL  = y(:,N+23*M:N+24*M-1);
model.TBGL    = y(:,N+24*M:N+25*M-1);
model.TTRL    = y(:,N+25*M:N+26*M-1);
model.ALBL    = y(:,N+26*M:N+27*M-1);
model.T4LT    = y(:,N+27*M:N+28*M-1);
model.T3LT    = y(:,N+28*M:N+29*M-1);
model.XL      = y(:,N+29*M:N+30*M-1);
model.XTTRL   = y(:,N+30*M:N+31*M-1);
model.XTBGL   = y(:,N+31*M:N+32*M-1);


%Free T4 in Liver and RB blood (Fig. S4A)
figure(401)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT4 in tissue blood with ALB only")

%T4ALB in Liver and RB blood (Fig. S4B)
figure(402)
hold on
plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4ALB in tissue blood with ALB only")

%Free T3 in Liver and RB blood (Fig. S4C)
figure(403)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT3 in tissue blood with ALB only")

%T3ALB in Liver and RB blood (Fig. S4D)
figure(404)
hold on
plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3ALB in tissue blood with ALB only")

%% ------------------------- Plotting Results (Fig. S5) --------------------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Setting TBG to zero 
param.TBGtot = 0;   

% Update the initial THBPs concentrations accordingly
init.TBGB    = param.TBGtot/param.VB;

y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600];  %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode15s('TH_PBK_Spatial_ODE',tspan, y0, options, param); 

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.XB      = y(:,23);
model.XTTRB   = y(:,24);
model.XT      = y(:,25);
model.XTTRT   = y(:,26);
model.XTBGB   = y(:,27);
model.XTBGT   = y(:,28);

N = 29;  %Starting index of repeating variables
M = param.M;
model.fT4RB   = y(:,N:N+M-1);
model.T4TBGRB = y(:,N+M:N+2*M-1);
model.T4TTRRB = y(:,N+2*M:N+3*M-1);
model.T4ALBRB = y(:,N+3*M:N+4*M-1);
model.fT3RB   = y(:,N+4*M:N+5*M-1);
model.T3TBGRB = y(:,N+5*M:N+6*M-1);
model.T3TTRRB = y(:,N+6*M:N+7*M-1);
model.T3ALBRB = y(:,N+7*M:N+8*M-1);
model.TBGRB   = y(:,N+8*M:N+9*M-1);
model.TTRRB   = y(:,N+9*M:N+10*M-1);
model.ALBRB   = y(:,N+10*M:N+11*M-1);
model.T4RBT   = y(:,N+11*M:N+12*M-1);
model.T3RBT   = y(:,N+12*M:N+13*M-1);
model.XRB     = y(:,N+13*M:N+14*M-1);
model.XTTRRB  = y(:,N+14*M:N+15*M-1);
model.XTBGRB  = y(:,N+15*M:N+16*M-1);

model.fT4L    = y(:,N+16*M:N+17*M-1);
model.T4TBGL  = y(:,N+17*M:N+18*M-1);
model.T4TTRL  = y(:,N+18*M:N+19*M-1);
model.T4ALBL  = y(:,N+19*M:N+20*M-1);
model.fT3L    = y(:,N+20*M:N+21*M-1);
model.T3TBGL  = y(:,N+21*M:N+22*M-1);
model.T3TTRL  = y(:,N+22*M:N+23*M-1);
model.T3ALBL  = y(:,N+23*M:N+24*M-1);
model.TBGL    = y(:,N+24*M:N+25*M-1);
model.TTRL    = y(:,N+25*M:N+26*M-1);
model.ALBL    = y(:,N+26*M:N+27*M-1);
model.T4LT    = y(:,N+27*M:N+28*M-1);
model.T3LT    = y(:,N+28*M:N+29*M-1);
model.XL      = y(:,N+29*M:N+30*M-1);
model.XTTRL   = y(:,N+30*M:N+31*M-1);
model.XTBGL   = y(:,N+31*M:N+32*M-1);


%Free T4 in Liver and RB blood (Fig. S5A)
figure(501)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT4 in tissue blood without TBG")

%T4TTR in Liver and RB blood (Fig. S5B)
figure(502)
hold on
plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4TTR in tissue blood without TBG")

%T4ALB in Liver and RB blood (Fig. S5C)
figure(503)
hold on
plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4ALB in tissue blood without TBG")


%Total T4 in Liver and RB blood (Fig. S5D)
figure(504)
hold on
plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4L(end,:)+model.T4TBGL(end,:)+model.T4TTRL(end,:)+model.T4ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4RB(end,:)+model.T4TBGRB(end,:)+model.T4TTRRB(end,:)+model.T4ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("Total T4 in tissue blood without TBG")

%Free T3 in Liver and RB blood (Fig. S5E)
figure(505)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT3 in tissue blood without TBG")

%T3TTR in Liver and RB blood (Fig. S5F)
figure(506)
hold on
plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3TTR in tissue blood without TBG")

%T3ALB in Liver and RB blood (Fig. S5G)
figure(507)
hold on
plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3ALB in tissue blood without TBG")

%Total T3 in Liver and RB blood (Fig. S5H)
figure(508)
hold on
plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3L(end,:)+model.T3TBGL(end,:)+model.T3TTRL(end,:)+model.T3ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3RB(end,:)+model.T3TBGRB(end,:)+model.T3TTRRB(end,:)+model.T3ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("Total T3 in tissue blood without TBG")


%% ------------------------- Plotting Results (Fig. S6) --------------------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Setting TTR to zero 
param.TTRtot = 0;   

% Update the initial THBPs concentrations accordingly
init.TTRB    = param.TTRtot/param.VB;

y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600];  %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode15s('TH_PBK_Spatial_ODE',tspan, y0, options, param); 

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.XB      = y(:,23);
model.XTTRB   = y(:,24);
model.XT      = y(:,25);
model.XTTRT   = y(:,26);
model.XTBGB   = y(:,27);
model.XTBGT   = y(:,28);

N = 29;  %Starting index of repeating variables
M = param.M;
model.fT4RB   = y(:,N:N+M-1);
model.T4TBGRB = y(:,N+M:N+2*M-1);
model.T4TTRRB = y(:,N+2*M:N+3*M-1);
model.T4ALBRB = y(:,N+3*M:N+4*M-1);
model.fT3RB   = y(:,N+4*M:N+5*M-1);
model.T3TBGRB = y(:,N+5*M:N+6*M-1);
model.T3TTRRB = y(:,N+6*M:N+7*M-1);
model.T3ALBRB = y(:,N+7*M:N+8*M-1);
model.TBGRB   = y(:,N+8*M:N+9*M-1);
model.TTRRB   = y(:,N+9*M:N+10*M-1);
model.ALBRB   = y(:,N+10*M:N+11*M-1);
model.T4RBT   = y(:,N+11*M:N+12*M-1);
model.T3RBT   = y(:,N+12*M:N+13*M-1);
model.XRB     = y(:,N+13*M:N+14*M-1);
model.XTTRRB  = y(:,N+14*M:N+15*M-1);
model.XTBGRB  = y(:,N+15*M:N+16*M-1);

model.fT4L    = y(:,N+16*M:N+17*M-1);
model.T4TBGL  = y(:,N+17*M:N+18*M-1);
model.T4TTRL  = y(:,N+18*M:N+19*M-1);
model.T4ALBL  = y(:,N+19*M:N+20*M-1);
model.fT3L    = y(:,N+20*M:N+21*M-1);
model.T3TBGL  = y(:,N+21*M:N+22*M-1);
model.T3TTRL  = y(:,N+22*M:N+23*M-1);
model.T3ALBL  = y(:,N+23*M:N+24*M-1);
model.TBGL    = y(:,N+24*M:N+25*M-1);
model.TTRL    = y(:,N+25*M:N+26*M-1);
model.ALBL    = y(:,N+26*M:N+27*M-1);
model.T4LT    = y(:,N+27*M:N+28*M-1);
model.T3LT    = y(:,N+28*M:N+29*M-1);
model.XL      = y(:,N+29*M:N+30*M-1);
model.XTTRL   = y(:,N+30*M:N+31*M-1);
model.XTBGL   = y(:,N+31*M:N+32*M-1);


%Free T4 in Liver and RB blood (Fig. S6A)
figure(601)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT4 in tissue blood without TTR")

%T4TBG in Liver and RB blood (Fig. S6B)
figure(602)
hold on
plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4TBG in tissue blood without TTR")

%T4ALB in Liver and RB blood (Fig. S6C)
figure(603)
hold on
plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4ALB in tissue blood without TTR")


%Total T4 in Liver and RB blood (Fig. S6D)
figure(604)
hold on
plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4L(end,:)+model.T4TBGL(end,:)+model.T4TTRL(end,:)+model.T4ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4RB(end,:)+model.T4TBGRB(end,:)+model.T4TTRRB(end,:)+model.T4ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("Total T4 in tissue blood without TTR")

%Free T3 in Liver and RB blood (Fig. S6E)
figure(605)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT3 in tissue blood without TTR")

%T3TBG in Liver and RB blood (Fig. S6F)
figure(606)
hold on
plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3TBG in tissue blood without TTR")

%T3ALB in Liver and RB blood (Fig. S6G)
figure(607)
hold on
plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3ALB in tissue blood without TTR")

%Total T3 in Liver and RB blood (Fig. S6H)
figure(608)
hold on
plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3L(end,:)+model.T3TBGL(end,:)+model.T3TTRL(end,:)+model.T3ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3RB(end,:)+model.T3TBGRB(end,:)+model.T3TTRRB(end,:)+model.T3ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("Total T3 in tissue blood without TTR")


%% ------------------------- Plotting Results (Fig. S7) --------------------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Setting ALB to zero 
param.ALBtot = 0;   

% Update the initial THBPs concentrations accordingly
init.ALBB    = param.ALBtot/param.VB;

y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600];  %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode15s('TH_PBK_Spatial_ODE',tspan, y0, options, param); 

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.XB      = y(:,23);
model.XTTRB   = y(:,24);
model.XT      = y(:,25);
model.XTTRT   = y(:,26);
model.XTBGB   = y(:,27);
model.XTBGT   = y(:,28);

N = 29;  %Starting index of repeating variables
M = param.M;
model.fT4RB   = y(:,N:N+M-1);
model.T4TBGRB = y(:,N+M:N+2*M-1);
model.T4TTRRB = y(:,N+2*M:N+3*M-1);
model.T4ALBRB = y(:,N+3*M:N+4*M-1);
model.fT3RB   = y(:,N+4*M:N+5*M-1);
model.T3TBGRB = y(:,N+5*M:N+6*M-1);
model.T3TTRRB = y(:,N+6*M:N+7*M-1);
model.T3ALBRB = y(:,N+7*M:N+8*M-1);
model.TBGRB   = y(:,N+8*M:N+9*M-1);
model.TTRRB   = y(:,N+9*M:N+10*M-1);
model.ALBRB   = y(:,N+10*M:N+11*M-1);
model.T4RBT   = y(:,N+11*M:N+12*M-1);
model.T3RBT   = y(:,N+12*M:N+13*M-1);
model.XRB     = y(:,N+13*M:N+14*M-1);
model.XTTRRB  = y(:,N+14*M:N+15*M-1);
model.XTBGRB  = y(:,N+15*M:N+16*M-1);

model.fT4L    = y(:,N+16*M:N+17*M-1);
model.T4TBGL  = y(:,N+17*M:N+18*M-1);
model.T4TTRL  = y(:,N+18*M:N+19*M-1);
model.T4ALBL  = y(:,N+19*M:N+20*M-1);
model.fT3L    = y(:,N+20*M:N+21*M-1);
model.T3TBGL  = y(:,N+21*M:N+22*M-1);
model.T3TTRL  = y(:,N+22*M:N+23*M-1);
model.T3ALBL  = y(:,N+23*M:N+24*M-1);
model.TBGL    = y(:,N+24*M:N+25*M-1);
model.TTRL    = y(:,N+25*M:N+26*M-1);
model.ALBL    = y(:,N+26*M:N+27*M-1);
model.T4LT    = y(:,N+27*M:N+28*M-1);
model.T3LT    = y(:,N+28*M:N+29*M-1);
model.XL      = y(:,N+29*M:N+30*M-1);
model.XTTRL   = y(:,N+30*M:N+31*M-1);
model.XTBGL   = y(:,N+31*M:N+32*M-1);


%Free T4 in Liver and RB blood (Fig. S7A)
figure(701)
hold on
plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT4 in tissue blood without ALB")

%T4TBG in Liver and RB blood (Fig. S7B)
figure(702)
hold on
plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4TBG in tissue blood without ALB")

%T4TTR in Liver and RB blood (Fig. S7C)
figure(703)
hold on
plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T4TTR in tissue blood without ALB")

%Total T4 in Liver and RB blood (Fig. S7D)
figure(704)
hold on
plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4L(end,:)+model.T4TBGL(end,:)+model.T4TTRL(end,:)+model.T4ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4RB(end,:)+model.T4TBGRB(end,:)+model.T4TTRRB(end,:)+model.T4ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("Total T4 in tissue blood without ALB")

%Free T3 in Liver and RB blood (Fig. S7E)
figure(705)
hold on
plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("fT3 in tissue blood without ALB")

%T3TBG in Liver and RB blood (Fig. S7F)
figure(706)
hold on
plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3TBG in tissue blood without ALB")

%T3TTR in Liver and RB blood (Fig. S7G)
figure(707)
hold on
plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("T3TTR in tissue blood without ALB")

%Total T3 in Liver and RB blood (Fig. S7H)
figure(708)
hold on
plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3L(end,:)+model.T3TBGL(end,:)+model.T3TTRL(end,:)+model.T3ALBL(end,:)], 'Color', [0.85,0.33,0.1], 'LineWidth', 1, 'LineStyle', '-')
plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3RB(end,:)+model.T3TBGRB(end,:)+model.T3TTRRB(end,:)+model.T3ALBRB(end,:)], 'Color', [0.47,0.67,0.19], 'LineWidth', 1, 'LineStyle', '-')
xlabel('Tissue Segment')
ylabel('Plasma Concentration (pM)')
legend("Liver","RB")
title("Total T3 in tissue blood without ALB")


%% ------------------------- Plotting Results (Fig. 8 and 9) ---------------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Clamp intial values of variables in Body Blood to basal steady-state values with corresponding ODEs set to zero as in TH_PBK_Spatial_ODE_Fig8.m
init.fT4B    = freeT4B;
init.T4TBGB  = T4TBGB;
init.T4TTRB  = T4TTRB;
init.T4ALBB  = T4ALBB;
init.fT3B    = freeT3B;
init.T3TBGB  = T3TBGB;
init.T3TTRB  = T3TTRB;
init.T3ALBB  = T3ALBB;
init.TBGB    = TBGB;
init.TTRB    = TTRB;
init.ALBB    = ALBB;


for caseIndex = 1:2 %1:5  % Iterate through 5 cases

    % Reset parameter to default
    param = default_param;

    % Use a switch statement for each case
    switch caseIndex
        case 1
            a2 = 0;
            param.k24 = 1.7084E-6/param.fuT4RBT*a2;
            param.k32 = 1.7084E-6/param.fuT4RBT*(1-a2);
            param.QRB  = 4 * QRBC*param.QC; %Plasma flow rate to RB
            color = [0.64 0.08 0.18];
            line_style = '-';
            display_name = 'k24=0, QRB=4xdefault';
        
        case 2
            a2 = 0;
            param.k24 = 1.7084E-6/param.fuT4RBT*a2;
            param.k32 = 1.7084E-6/param.fuT4RBT*(1-a2);
            param.QRB  = 0.25 * QRBC*param.QC; %Plasma flow rate to RB
            color = [0.64 0.08 0.18];
            line_style = '--';
            display_name = 'k24=0, QRB=0.25xdefault';
       case 3
            param.k33 = 0;
            param.QRB  = 4 * QRBC*param.QC; %Plasma flow rate to RB
            color = [0.30 0.75 0.93];
            line_style = '-';
            display_name = 'k33=0, QRB=4xdefault';
        
        case 4
            param.k33 = 0;
            param.QRB  = 0.25 * QRBC*param.QC; %Plasma flow rate to RB
            color = [0.30 0.75 0.93];
            line_style = '--';
            display_name = 'k33=0, QRB=0.25xdefault';
        
        case 5
            color = [0.00 0.00 0.00];
            line_style = '-';
            display_name = 'Default';
    end
    
    y0 = cell2mat(struct2cell(init));
    tspan = [0:10000:1000*24*3600];  %Running for 1000 days
    options = odeset('RelTol',1e-8,'AbsTol',1e-6);
    [t,y] = ode15s('TH_PBK_Spatial_ODE_Fig8',tspan, y0, options, param); 
    
    %Simulation results
    model.fT4B    = y(:,1);
    model.T4TBGB  = y(:,2);
    model.T4TTRB  = y(:,3);
    model.T4ALBB  = y(:,4);
    
    model.fT3B    = y(:,5);
    model.T3TBGB  = y(:,6);
    model.T3TTRB  = y(:,7);
    model.T3ALBB  = y(:,8);
    
    model.TBGB    = y(:,9);
    model.TTRB    = y(:,10);
    model.ALBB    = y(:,11);
    
    model.fT4T    = y(:,12);
    model.T4TBGT  = y(:,13);
    model.T4TTRT  = y(:,14);
    model.T4ALBT  = y(:,15);
    model.fT3T    = y(:,16);
    model.T3TBGT  = y(:,17);
    model.T3TTRT  = y(:,18);
    model.T3ALBT  = y(:,19);
    model.TBGT    = y(:,20);
    model.TTRT    = y(:,21);
    model.ALBT    = y(:,22);
    
    model.XB      = y(:,23);
    model.XTTRB   = y(:,24);
    model.XT      = y(:,25);
    model.XTTRT   = y(:,26);
    model.XTBGB   = y(:,27);
    model.XTBGT   = y(:,28);
    
    N = 29;  %Starting index of repeating variables
    M = param.M;
    model.fT4RB   = y(:,N:N+M-1);
    model.T4TBGRB = y(:,N+M:N+2*M-1);
    model.T4TTRRB = y(:,N+2*M:N+3*M-1);
    model.T4ALBRB = y(:,N+3*M:N+4*M-1);
    model.fT3RB   = y(:,N+4*M:N+5*M-1);
    model.T3TBGRB = y(:,N+5*M:N+6*M-1);
    model.T3TTRRB = y(:,N+6*M:N+7*M-1);
    model.T3ALBRB = y(:,N+7*M:N+8*M-1);
    model.TBGRB   = y(:,N+8*M:N+9*M-1);
    model.TTRRB   = y(:,N+9*M:N+10*M-1);
    model.ALBRB   = y(:,N+10*M:N+11*M-1);
    model.T4RBT   = y(:,N+11*M:N+12*M-1);
    model.T3RBT   = y(:,N+12*M:N+13*M-1);
    model.XRB     = y(:,N+13*M:N+14*M-1);
    model.XTTRRB  = y(:,N+14*M:N+15*M-1);
    model.XTBGRB  = y(:,N+15*M:N+16*M-1);
    
    model.fT4L    = y(:,N+16*M:N+17*M-1);
    model.T4TBGL  = y(:,N+17*M:N+18*M-1);
    model.T4TTRL  = y(:,N+18*M:N+19*M-1);
    model.T4ALBL  = y(:,N+19*M:N+20*M-1);
    model.fT3L    = y(:,N+20*M:N+21*M-1);
    model.T3TBGL  = y(:,N+21*M:N+22*M-1);
    model.T3TTRL  = y(:,N+22*M:N+23*M-1);
    model.T3ALBL  = y(:,N+23*M:N+24*M-1);
    model.TBGL    = y(:,N+24*M:N+25*M-1);
    model.TTRL    = y(:,N+25*M:N+26*M-1);
    model.ALBL    = y(:,N+26*M:N+27*M-1);
    model.T4LT    = y(:,N+27*M:N+28*M-1);
    model.T3LT    = y(:,N+28*M:N+29*M-1);
    model.XL      = y(:,N+29*M:N+30*M-1);
    model.XTTRL   = y(:,N+30*M:N+31*M-1);
    model.XTBGL   = y(:,N+31*M:N+32*M-1);
    

    %Free T4 in RB blood (Fig. 8A)
    figure(81)
    hold on
    plot(0:1:param.M,[model.fT4B(end),model.fT4RB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood fT4")
    legend show;
    
    %T4TBG in RB blood (Fig. 8B)
    figure(82)
    hold on
    plot(0:1:param.M,[model.T4TBGB(end),model.T4TBGRB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood T4TBG")
    legend show;
    
    %T4TTR in RB blood (Fig. 8C)
    figure(83)
    hold on
    plot(0:1:param.M,[model.T4TTRB(end),model.T4TTRRB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood T4TTR")
    legend show;
    
    %T4ALB in RB blood (Fig. 8D)
    figure(84)
    hold on
    plot(0:1:param.M,[model.T4ALBB(end),model.T4ALBRB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood T4ALB")
    legend show;
    
    %Total T4 in RB blood (Fig. 8E)
    figure(85)
    hold on
    plot(0:1:param.M,[model.fT4B(end) + model.T4TBGB(end) + model.T4TTRB(end) + model.T4ALBB(end), model.fT4RB(end,:) + model.T4TBGRB(end,:) + model.T4TTRRB(end,:) + model.T4ALBRB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood Total T4")
    legend show;
    

    %Free T3 in RB blood (Fig. 8F)
    figure(86)
    hold on
    plot(0:1:param.M,[model.fT3B(end),model.fT3RB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood fT3")
    legend show;
    
    %T3TBG in RB blood (Fig. 8G)
    figure(87)
    hold on
    plot(0:1:param.M,[model.T3TBGB(end),model.T3TBGRB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood T3TBG")
    legend show;
    
    %T3TTR in RB blood (Fig. 8H)
    figure(88)
    hold on
    plot(0:1:param.M,[model.T3TTRB(end),model.T3TTRRB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood T3TTR")
    legend show;
    
    %T3ALB in RB blood (Fig. 8I)
    figure(89)
    hold on
    plot(0:1:param.M,[model.T3ALBB(end),model.T3ALBRB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood T3ALB")
    legend show;
    
    %Total T3 in RB blood (Fig. 8J)
    figure(90)
    hold on
    plot(0:1:param.M,[model.fT3B(end) + model.T3TBGB(end) + model.T3TTRB(end) + model.T3ALBB(end), model.fT3RB(end,:) + model.T3TBGRB(end,:) + model.T3TTRRB(end,:) + model.T3ALBRB(end,:)], 'Color', color, 'LineWidth', 1, 'LineStyle', line_style, 'DisplayName', display_name)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration (pM)')
    title("RB blood Total T3")
    legend show;


    %Absolute Drop of T4 in RB blood (Fig. 9A, 9C, 9E, and 9G)
    figure(90+caseIndex*2-1)
    hold on
    plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBRB(end,1) + model.T4TTRB(end)-model.T4TTRRB(end,1) + model.T4TBGB(end)-model.T4TBGRB(end,1) + model.fT4B(end)-model.fT4RB(end,1), - model.T4TBGRB(end,2:param.M) + model.T4TBGRB(end,1:param.M-1) - model.T4ALBRB(end,2:param.M) + model.T4ALBRB(end,1:param.M-1) - model.T4TTRRB(end,2:param.M) + model.T4TTRRB(end,1:param.M-1) - model.fT4RB(end,2:param.M) + model.fT4RB(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
    plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBRB(end,1), - model.T4ALBRB(end,2:param.M) + model.T4ALBRB(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
    plot([1:1:param.M], [model.T4TTRB(end)-model.T4TTRRB(end,1), - model.T4TTRRB(end,2:param.M) + model.T4TTRRB(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
    plot([1:1:param.M], [model.T4TBGB(end)-model.T4TBGRB(end,1), - model.T4TBGRB(end,2:param.M) + model.T4TBGRB(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
    plot([1:1:param.M], [model.fT4B(end)-model.fT4RB(end,1), - model.fT4RB(end,2:param.M) + model.fT4RB(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration Drop (pM)')
    legend("Total T4", "T4ALB", "T4TTR", "T4TBG", "fT4")
    title(sprintf('RB blood T4: %s', display_name))
    
    %Absolute Drop of T3 in RB blood (Fig. 9B, 9D, 9F, and 9H)
    figure(91+caseIndex*2-1)
    hold on
    plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBRB(end,1) + model.T3TTRB(end)-model.T3TTRRB(end,1) + model.T3TBGB(end)-model.T3TBGRB(end,1) + model.fT3B(end)-model.fT3RB(end,1), - model.T3TBGRB(end,2:param.M) + model.T3TBGRB(end,1:param.M-1) - model.T3ALBRB(end,2:param.M) + model.T3ALBRB(end,1:param.M-1) - model.T3TTRRB(end,2:param.M) + model.T3TTRRB(end,1:param.M-1) - model.fT3RB(end,2:param.M) + model.fT3RB(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
    plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBRB(end,1), - model.T3ALBRB(end,2:param.M) + model.T3ALBRB(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
    plot([1:1:param.M], [model.T3TTRB(end)-model.T3TTRRB(end,1), - model.T3TTRRB(end,2:param.M) + model.T3TTRRB(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
    plot([1:1:param.M], [model.T3TBGB(end)-model.T3TBGRB(end,1), - model.T3TBGRB(end,2:param.M) + model.T3TBGRB(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
    plot([1:1:param.M], [model.fT3B(end)-model.fT3RB(end,1), - model.fT3RB(end,2:param.M) + model.fT3RB(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
    xlabel('Tissue Segment')
    ylabel('Plasma Concentration Drop (pM)')
    legend("Total T3", "T3ALB", "T3TTR", "T3TBG", "fT3")
    title(sprintf('RB blood T3: %s', display_name))

end


%% ------------------------- Plotting Results (Fig. 10-12, S8-S11) ---------------------- %%
% Reset parameter and initial condition to default
param = default_param;
init = init_default;

% Clamp intial values of variables in Body Blood to basal steady-state values with corresponding ODEs set to zero as in TH_PBK_Spatial_ODE_gradient.m
init.fT4B    = freeT4B;
init.T4TBGB  = T4TBGB;
init.T4TTRB  = T4TTRB;
init.T4ALBB  = T4ALBB;
init.fT3B    = freeT3B;
init.T3TBGB  = T3TBGB;
init.T3TTRB  = T3TTRB;
init.T3ALBB  = T3ALBB;
init.TBGB    = TBGB;
init.TTRB    = TTRB;
init.ALBB    = ALBB;

% List of parameters with gradient
arrangement_of_values = 1:1:param.M; 
parameters = ["k25" "k30" "k27" "k31" "k26" "k34" "k35" "k21" "k28" "k23" "k29" "k24" "k32" "k33"];
parameter_indices = 1:length(parameters);
parameter_dictionary = containers.Map(parameter_indices, parameters);

% Running model with one parameter gradient at a time
for param_index = 0:14 % when param_index=0, it runs the default case (no gradient); param_index=1-7, it runs the gradient in Liver; param_index=8-14, it runs the gradient in RB; 

    param_index

    % Increasing or decreasing gradient when arrangement=1 and 2 respectively.
    for arrangement = 1:2

        tic
        
        param.k25 = arrayfun(@(current) get_gradient_value(param_index, 1, default_param.k25, current, param.M), arrangement_of_values); % influx_T4_liver_gradients 
        
        param.k30 = arrayfun(@(current) get_gradient_value(param_index, 2, default_param.k30, current, param.M), arrangement_of_values); % efflux_T4_liver_gradients
        
        param.k27 = arrayfun(@(current) get_gradient_value(param_index, 3, default_param.k27, current, param.M), arrangement_of_values); % influx_T3_liver_gradients
        
        param.k31 = arrayfun(@(current) get_gradient_value(param_index, 4, default_param.k31, current, param.M), arrangement_of_values); % efflux_T3_liver_gradient 
        
        param.k26 = arrayfun(@(current) get_gradient_value(param_index, 5, default_param.k26, current, param.M), arrangement_of_values); % conversion_liver_gradients
        
        param.k34 = arrayfun(@(current) get_gradient_value(param_index, 6, default_param.k34, current, param.M), arrangement_of_values); % degradation_T4_liver_gradients
        
        param.k35 = arrayfun(@(current) get_gradient_value(param_index, 7, default_param.k35, current, param.M), arrangement_of_values); % degradation_T3_liver_gradients

        param.k21 = arrayfun(@(current) get_gradient_value(param_index, 8, default_param.k21, current, param.M), arrangement_of_values); % influx_T4_RB_gradients

        param.k28 = arrayfun(@(current) get_gradient_value(param_index, 9, default_param.k28, current, param.M), arrangement_of_values); % efflux_T4_RB_gradients

        param.k23 = arrayfun(@(current) get_gradient_value(param_index, 10, default_param.k23, current, param.M), arrangement_of_values); % influx_T3_RB_gradients

        param.k29 = arrayfun(@(current) get_gradient_value(param_index, 11, default_param.k29, current, param.M), arrangement_of_values); % efflux_T3_RB_gradients

        param.k24 = arrayfun(@(current) get_gradient_value(param_index, 12, default_param.k24, current, param.M), arrangement_of_values); % conversion_RB_gradients

        param.k32 = arrayfun(@(current) get_gradient_value(param_index, 13, default_param.k32, current, param.M), arrangement_of_values); % degradation_T4_RB_gradients
        
        param.k33 = arrayfun(@(current) get_gradient_value(param_index, 14, default_param.k33, current, param.M), arrangement_of_values); % degradation_T3_RB_gradients

        % ------------------------- Run To Steady-State ----------------------- %
        y0 = cell2mat(struct2cell(init));
        tspan = [0:10000:1000*24*3600];  %Running for 1000 days
        options = odeset('RelTol',1e-8,'AbsTol',1e-6);
        [t,y] = ode15s('TH_PBK_Spatial_ODE_gradient',tspan, y0, options, param); 
        
        %Simulation results
        model.fT4B    = y(:,1);
        model.T4TBGB  = y(:,2);
        model.T4TTRB  = y(:,3);
        model.T4ALBB  = y(:,4);
        
        model.fT3B    = y(:,5);
        model.T3TBGB  = y(:,6);
        model.T3TTRB  = y(:,7);
        model.T3ALBB  = y(:,8);
        
        model.TBGB    = y(:,9);
        model.TTRB    = y(:,10);
        model.ALBB    = y(:,11);
        
        model.fT4T    = y(:,12);
        model.T4TBGT  = y(:,13);
        model.T4TTRT  = y(:,14);
        model.T4ALBT  = y(:,15);
        model.fT3T    = y(:,16);
        model.T3TBGT  = y(:,17);
        model.T3TTRT  = y(:,18);
        model.T3ALBT  = y(:,19);
        model.TBGT    = y(:,20);
        model.TTRT    = y(:,21);
        model.ALBT    = y(:,22);
        
        model.XB      = y(:,23);
        model.XTTRB   = y(:,24);
        model.XT      = y(:,25);
        model.XTTRT   = y(:,26);
        model.XTBGB   = y(:,27);
        model.XTBGT   = y(:,28);
        
        N = 29;  %Starting index of repeating variables
        M = param.M;
        model.fT4RB   = y(:,N:N+M-1);
        model.T4TBGRB = y(:,N+M:N+2*M-1);
        model.T4TTRRB = y(:,N+2*M:N+3*M-1);
        model.T4ALBRB = y(:,N+3*M:N+4*M-1);
        model.fT3RB   = y(:,N+4*M:N+5*M-1);
        model.T3TBGRB = y(:,N+5*M:N+6*M-1);
        model.T3TTRRB = y(:,N+6*M:N+7*M-1);
        model.T3ALBRB = y(:,N+7*M:N+8*M-1);
        model.TBGRB   = y(:,N+8*M:N+9*M-1);
        model.TTRRB   = y(:,N+9*M:N+10*M-1);
        model.ALBRB   = y(:,N+10*M:N+11*M-1);
        model.T4RBT   = y(:,N+11*M:N+12*M-1);
        model.T3RBT   = y(:,N+12*M:N+13*M-1);
        model.XRB     = y(:,N+13*M:N+14*M-1);
        model.XTTRRB  = y(:,N+14*M:N+15*M-1);
        model.XTBGRB  = y(:,N+15*M:N+16*M-1);
        
        model.fT4L    = y(:,N+16*M:N+17*M-1);
        model.T4TBGL  = y(:,N+17*M:N+18*M-1);
        model.T4TTRL  = y(:,N+18*M:N+19*M-1);
        model.T4ALBL  = y(:,N+19*M:N+20*M-1);
        model.fT3L    = y(:,N+20*M:N+21*M-1);
        model.T3TBGL  = y(:,N+21*M:N+22*M-1);
        model.T3TTRL  = y(:,N+22*M:N+23*M-1);
        model.T3ALBL  = y(:,N+23*M:N+24*M-1);
        model.TBGL    = y(:,N+24*M:N+25*M-1);
        model.TTRL    = y(:,N+25*M:N+26*M-1);
        model.ALBL    = y(:,N+26*M:N+27*M-1);
        model.T4LT    = y(:,N+27*M:N+28*M-1);
        model.T3LT    = y(:,N+28*M:N+29*M-1);
        model.XL      = y(:,N+29*M:N+30*M-1);
        model.XTTRL   = y(:,N+30*M:N+31*M-1);
        model.XTBGL   = y(:,N+31*M:N+32*M-1);
        
        
  % ------------------------- Plotting Results ---------------------------------------------------- %
        
    % -------------------Parameter Gradient Plot (Panels A in Fig.10-12,S8-S11)-----%
        figure(1000 + 100 * (param_index - 1))
        hold on
        
        %k25 (Fig. 10A)
        if param_index == 1
             if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k25, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k25, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k25, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (L/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
             end
        
        %k30 (Fig. 11A)
        elseif param_index == 2
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k30, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k30, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k30, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (L/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k27 (Fig. S8A)
        elseif param_index == 3
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k27, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k27, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k27, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (L/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k31 (Fig. S9A)
        elseif param_index == 4
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k31, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k31, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k31, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (L/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k26 (Fig. 12A)
        elseif param_index == 5
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k26, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k26, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k26, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (1/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k34 (Fig. S10A)
        elseif param_index == 6
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k34, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k34, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k34, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (1/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k35 (Fig. S11A)
        elseif param_index == 7
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k35, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k35, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k35, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (1/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k21 
        elseif param_index == 8
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k21, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k21, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k21, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (L/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k28 
        elseif param_index == 9
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k28, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k28, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k28, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (L/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k23 
        elseif param_index == 10
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k23, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k23, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k23, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (L/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k29 
        elseif param_index == 11
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k29, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k29, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k29, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (L/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k24 
        elseif param_index == 12
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k24, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k24, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k24, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (1/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k32 
        elseif param_index == 13
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k32, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k32, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k32, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (1/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %k33 
        elseif param_index == 14
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k33, 1:1:param.M), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
                plot(1:1:param.M, param.k33, 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
            else
                plot(1:1:param.M, param.k33, 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
                xlabel('Tissue Segment')
                ylabel(sprintf('%s (1/S)', parameter_dictionary(param_index)))
                legend show
                title(sprintf('%s', parameter_dictionary(param_index)))
            end
        
        %default (no gradient) 
        else
            default_model = model;
            break
        end


 % --------------------Liver----------------------------------------------------------------------------------%       
   if param_index >=1 && param_index <= 7
   % -------------------Free T4 in Liver blood (Panels B in Fig.10-12,S10) -----%
        figure(1000 + 100 * (param_index - 1) + 1)
        hold on
        if arrangement == 1
            plot([0:1:param.M],[default_model.fT4B(end),default_model.fT4RB(end,:)], 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)], 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Plasma Concentration (pM)')
            legend show
            title("Liver blood fT4")
        end

       
   % -------------------Total T4 in Liver blood (Panels C in Fig.10-12,S10) -----%
        figure(figure(1000 + 100 * (param_index - 1) + 2))
        hold on
        if arrangement == 1
            plot([0:1:param.M],[default_model.fT4B(end)+default_model.T4TBGB(end)+default_model.T4TTRB(end)+default_model.T4ALBB(end), default_model.fT4L(end,:)+default_model.T4TBGL(end,:)+default_model.T4TTRL(end,:)+default_model.T4ALBL(end,:)], 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4L(end,:)+model.T4TBGL(end,:)+model.T4TTRL(end,:)+model.T4ALBL(end,:)], 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4L(end,:)+model.T4TBGL(end,:)+model.T4TTRL(end,:)+model.T4ALBL(end,:)], 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Plasma Concentration (pM)')
            legend show
            title("Liver blood Total T4")
        end

        
    % -------------------Free T3 in Liver blood (Panels G in Fig. 10-12,S10, panels D in Fig. S8, S9, and S11) -----%
        figure(figure(1000 + 100 * (param_index - 1) + 6))
        hold on
        if arrangement == 1
            plot([0:1:param.M],[default_model.fT3B(end),default_model.fT3L(end,:)], 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)], 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Plasma Concentration (pM)')
            legend show
            title("Liver blood fT3")
        end

        
    % -------------------Total T3 in Liver blood (Panels H in Fig. 10-12 and S10, panels E in Fig. S8, S9, and S11) -----%
        figure(figure(1000 + 100 * (param_index - 1) + 7))
        % hold on
        % plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3L(end,:)+model.T3TBGL(end,:)+model.T3TTRL(end,:)+model.T3ALBL(end,:)])
        % plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3RB(end,:)+model.T3TBGRB(end,:)+model.T3TTRRB(end,:)+model.T3ALBRB(end,:)])
         
        hold on
        if arrangement == 1
            plot([0:1:param.M],[default_model.fT3B(end)+default_model.T3TBGB(end)+default_model.T3TTRB(end)+default_model.T3ALBB(end), default_model.fT3L(end,:)+default_model.T3TBGL(end,:)+default_model.T3TTRL(end,:)+default_model.T3ALBL(end,:)], 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3L(end,:)+model.T3TBGL(end,:)+model.T3TTRL(end,:)+model.T3ALBL(end,:)], 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3L(end,:)+model.T3TBGL(end,:)+model.T3TTRL(end,:)+model.T3ALBL(end,:)], 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Plasma Concentration (pM)')
            legend show
            title("Liver blood Total T3")
        end

        
    % -------------------T4 in Liver tissue (Panels D in Fig. 10-12 and S10) -----%
        figure(figure(1000 + 100 * (param_index - 1) + 3))
        % hold on
        % plot([1:1:param.M],model.T4LT(end,:))
        % plot([1:1:param.M],model.T4RBT(end,:))

        hold on
        if arrangement == 1
            plot([1:1:param.M],default_model.T4LT(end,:), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([1:1:param.M],model.T4LT(end,:), 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([1:1:param.M],model.T4LT(end,:), 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Concentration (pM)')
            legend show
            title("Liver tissue T4")
        end

        
    % -------------------T3 in Liver tissue (Panels I in Fig. 10-12 and S10, panels F in Fig. S8, S9, and S11) -----%
        figure(figure(1000 + 100 * (param_index - 1) + 8))
        % hold on
        % plot([1:1:param.M],model.T3LT(end,:))
        % plot([1:1:param.M],model.T3RBT(end,:))

        hold on
        if arrangement == 1
            plot([1:1:param.M],default_model.T3LT(end,:), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([1:1:param.M],model.T3LT(end,:), 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([1:1:param.M],model.T3LT(end,:), 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Concentration (pM)')
            legend show
            title("Liver tissue T3")
        end
        

        if arrangement == 1
            arrow = '↑';
        else
            arrow = '↓';
        end
    % -------------------Absolute drop of T4 in Liver blood (Panels E and F in Fig. 10-12 and S10) -----%
        figure(figure(1000 + 100 * (param_index - 1) + 4 + arrangement - 1))
        hold on
        plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBL(end,1) + model.T4TTRB(end)-model.T4TTRL(end,1) + model.T4TBGB(end)-model.T4TBGL(end,1) + model.fT4B(end)-model.fT4L(end,1), - model.T4TBGL(end,2:param.M) + model.T4TBGL(end,1:param.M-1) - model.T4ALBL(end,2:param.M) + model.T4ALBL(end,1:param.M-1) - model.T4TTRL(end,2:param.M) + model.T4TTRL(end,1:param.M-1) - model.fT4L(end,2:param.M) + model.fT4L(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBL(end,1), - model.T4ALBL(end,2:param.M) + model.T4ALBL(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T4TTRB(end)-model.T4TTRL(end,1), - model.T4TTRL(end,2:param.M) + model.T4TTRL(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T4TBGB(end)-model.T4TBGL(end,1), - model.T4TBGL(end,2:param.M) + model.T4TBGL(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.fT4B(end)-model.fT4L(end,1), - model.fT4L(end,2:param.M) + model.fT4L(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
        xlabel('Tissue Segment')
        ylabel('Plasma Concentration Drop (pM)')
        legend("Total T4", "T4ALB", "T4TTR", "T4TBG", "fT4")
        title(sprintf("Liver blood T4: %s %s", parameter_dictionary(param_index), arrow))
        
     % -------------------Absolute drop of T3 in Liver blood (Panels J and K in Fig. 10-12 and S10, panels B and C in Fig. S8, S9, and S11) -----%
        figure(1000 + 100 * (param_index - 1) + 9 + arrangement - 1)
        hold on
        plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBL(end,1) + model.T3TTRB(end)-model.T3TTRL(end,1) + model.T3TBGB(end)-model.T3TBGL(end,1) + model.fT3B(end)-model.fT3L(end,1), - model.T3TBGL(end,2:param.M) + model.T3TBGL(end,1:param.M-1) - model.T3ALBL(end,2:param.M) + model.T3ALBL(end,1:param.M-1) - model.T3TTRL(end,2:param.M) + model.T3TTRL(end,1:param.M-1) - model.fT3L(end,2:param.M) + model.fT3L(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBL(end,1), - model.T3ALBL(end,2:param.M) + model.T3ALBL(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T3TTRB(end)-model.T3TTRL(end,1), - model.T3TTRL(end,2:param.M) + model.T3TTRL(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T3TBGB(end)-model.T3TBGL(end,1), - model.T3TBGL(end,2:param.M) + model.T3TBGL(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.fT3B(end)-model.fT3L(end,1), - model.fT3L(end,2:param.M) + model.fT3L(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
        xlabel('Tissue Segment')
        ylabel('Plasma Concentration Drop (pM)')
        legend("Total T3", "T3ALB", "T3TTR", "T3TBG", "fT3")
        title(sprintf("Liver blood T3: %s %s", parameter_dictionary(param_index), arrow))
  


  % --------------------RB----------------------------------------------------------------------------------%       
    elseif param_index >= 8 && param_index <= 14
    
     % -------------------Free T4 in RB blood-----%
        figure(1000 + 100 * (param_index - 1) + 1)
        hold on
        if arrangement == 1
            plot([0:1:param.M],[default_model.fT4B(end),default_model.fT4RB(end,:)], 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)], 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Plasma Concentration (pM)')
            legend show
            title(sprintf('RB blood fT4: %s', parameter_dictionary(param_index)))
        end

     % -------------------Total T4 in RB blood -----%
        figure(figure(1000 + 100 * (param_index - 1) + 2))
        hold on
        if arrangement == 1
            plot([0:1:param.M],[default_model.fT4B(end)+default_model.T4TBGB(end)+default_model.T4TTRB(end)+default_model.T4ALBB(end), default_model.fT4RB(end,:)+default_model.T4TBGRB(end,:)+default_model.T4TTRRB(end,:)+default_model.T4ALBRB(end,:)], 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4RB(end,:)+model.T4TBGRB(end,:)+model.T4TTRRB(end,:)+model.T4ALBRB(end,:)], 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4RB(end,:)+model.T4TBGRB(end,:)+model.T4TTRRB(end,:)+model.T4ALBRB(end,:)], 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Plasma Concentration (pM)')
            legend show
            title(sprintf('RB blood Total T4: %s', parameter_dictionary(param_index)))
        end

     % -------------------Free T3 in RB blood -----%
        figure(figure(1000 + 100 * (param_index - 1) + 6))
        hold on
        if arrangement == 1
            plot([0:1:param.M],[default_model.fT3B(end),default_model.fT3RB(end,:)], 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)], 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Plasma Concentration (pM)')
            legend show
            title(sprintf('RB blood fT3: %s', parameter_dictionary(param_index)))
        end
        
     % -------------------Total T3 in RB blood -----%
        figure(figure(1000 + 100 * (param_index - 1) + 7))
        hold on
        if arrangement == 1
            plot([0:1:param.M],[default_model.fT3B(end)+default_model.T3TBGB(end)+default_model.T3TTRB(end)+default_model.T3ALBB(end), default_model.fT3RB(end,:)+default_model.T3TBGRB(end,:)+default_model.T3TTRRB(end,:)+default_model.T3ALBRB(end,:)], 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3RB(end,:)+model.T3TBGRB(end,:)+model.T3TTRRB(end,:)+model.T3ALBRB(end,:)], 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3RB(end,:)+model.T3TBGRB(end,:)+model.T3TTRRB(end,:)+model.T3ALBRB(end,:)], 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Plasma Concentration (pM)')
            legend show
            title(sprintf('RB blood Total T3: %s', parameter_dictionary(param_index)))
        end

     % -------------------T4 in RB tissue -----%
        figure(figure(1000 + 100 * (param_index - 1) + 3))
        hold on
        if arrangement == 1
            plot([1:1:param.M],default_model.T4RBT(end,:), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([1:1:param.M],model.T4RBT(end,:), 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([1:1:param.M],model.T4RBT(end,:), 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Concentration (pM)')
            legend show
            title(sprintf('RB tissue T4: %s', parameter_dictionary(param_index)))
        end

     % -------------------T3 in RB tissue -----%
        figure(figure(1000 + 100 * (param_index - 1) + 8))
        hold on
        if arrangement == 1
            plot([1:1:param.M],default_model.T3RBT(end,:), 'Color', [0.50 0.50 0.50], 'DisplayName', sprintf("%s ↔", parameter_dictionary(param_index)))
            plot([1:1:param.M],model.T3RBT(end,:), 'Color', [0.64 0.08 0.18], 'DisplayName', sprintf("%s ↑", parameter_dictionary(param_index)))
        else
            plot([1:1:param.M],model.T3RBT(end,:), 'Color', [0.47 0.67 0.19], 'DisplayName', sprintf("%s ↓", parameter_dictionary(param_index)))
            xlabel('Tissue Segment')
            ylabel('Concentration (pM)')
            legend show
            title(sprintf('RB tissue T3: %s', parameter_dictionary(param_index)))
        end
        

        if arrangement == 1
            arrow = '↑';
        else
            arrow = '↓';
        end
     % -------------------Absolute drop of T4 in RB blood -----%
        figure(figure(1000 + 100 * (param_index - 1) + 4 + arrangement - 1))
        hold on
        plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBRB(end,1) + model.T4TTRB(end)-model.T4TTRRB(end,1) + model.T4TBGB(end)-model.T4TBGRB(end,1) + model.fT4B(end)-model.fT4RB(end,1), - model.T4TBGRB(end,2:param.M) + model.T4TBGRB(end,1:param.M-1) - model.T4ALBRB(end,2:param.M) + model.T4ALBRB(end,1:param.M-1) - model.T4TTRRB(end,2:param.M) + model.T4TTRRB(end,1:param.M-1) - model.fT4RB(end,2:param.M) + model.fT4RB(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBRB(end,1), - model.T4ALBRB(end,2:param.M) + model.T4ALBRB(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T4TTRB(end)-model.T4TTRRB(end,1), - model.T4TTRRB(end,2:param.M) + model.T4TTRRB(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T4TBGB(end)-model.T4TBGRB(end,1), - model.T4TBGRB(end,2:param.M) + model.T4TBGRB(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.fT4B(end)-model.fT4RB(end,1), - model.fT4RB(end,2:param.M) + model.fT4RB(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
        xlabel('Tissue Segment')
        ylabel('Plasma Concentration Drop (pM)')
        legend("Total T4", "T4ALB", "T4TTR", "T4TBG", "fT4")
        title(sprintf("RB blood T4: %s %s", parameter_dictionary(param_index), arrow))
        
     % -------------------Absolute drop of T3 in RB blood -----%
        figure(1000 + 100 * (param_index - 1) + 9 + arrangement - 1)
        hold on
        plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBRB(end,1) + model.T3TTRB(end)-model.T3TTRRB(end,1) + model.T3TBGB(end)-model.T3TBGRB(end,1) + model.fT3B(end)-model.fT3RB(end,1), - model.T3TBGRB(end,2:param.M) + model.T3TBGRB(end,1:param.M-1) - model.T3ALBRB(end,2:param.M) + model.T3ALBRB(end,1:param.M-1) - model.T3TTRRB(end,2:param.M) + model.T3TTRRB(end,1:param.M-1) - model.fT3RB(end,2:param.M) + model.fT3RB(end,1:param.M-1)], 'Color', [0.39,0.83,0.07], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBRB(end,1), - model.T3ALBRB(end,2:param.M) + model.T3ALBRB(end,1:param.M-1)], 'Color', [0,0.45,0.74], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T3TTRB(end)-model.T3TTRRB(end,1), - model.T3TTRRB(end,2:param.M) + model.T3TTRRB(end,1:param.M-1)], 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.T3TBGB(end)-model.T3TBGRB(end,1), - model.T3TBGRB(end,2:param.M) + model.T3TBGRB(end,1:param.M-1)], 'Color', [0.93,0.69,0.13], 'LineWidth', 1.5)
        plot([1:1:param.M], [model.fT3B(end)-model.fT3RB(end,1), - model.fT3RB(end,2:param.M) + model.fT3RB(end,1:param.M-1)], 'Color', [0.64,0.08,0.18], 'LineWidth', 1.5)
        xlabel('Tissue Segment')
        ylabel('Plasma Concentration Drop (pM)')
        legend("Total T3", "T3ALB", "T3TTR", "T3TBG", "fT3")
        title(sprintf("RB blood T3: %s %s", parameter_dictionary(param_index), arrow))
        
    end

        arrangement_of_values = flip(arrangement_of_values);

        toc

    end
end


%% ------------------------- Sensitivity Analysis (Fig. 2C and 2D) -----------------------%%

y0 = y0_default_steady_state; %Using the steady-state values from the "Run to Steady State" section above as the initial values

sensitivity_analysis_params = ["k1", "k2", "k3", "k4", "k5", "k6", "k7"...
    "k8", "k9", "k10", "k11", "k12", "TBGtot", "TTRtot", "ALBtot", "QT", "QL", "QRB"];

sensitivity_analysis_results = struct();

%Simulation results

N = 29;  %Starting index of repeating variables
M = param.M;

% T4 and T3 differnetial between arterial and venous blood concentrations
ss_freeT4B = y0(end,1);
ss_T4TBGB = y0(end,2);
ss_T4TTRB = y0(end,3);
ss_T4ALBB = y0(end,4);

ss_freeT3B = y0(end,5);
ss_T3TBGB = y0(end,6);
ss_T3TTRB = y0(end,7);
ss_T3ALBB = y0(end,8);

% Thyroid blood Compartment
ss_freeT4T = y0(end,12);
ss_T4TBGT = y0(end,13);
ss_T4TTRT = y0(end,14);
ss_T4ALBT = y0(end,15);

ss_freeT3T = y0(end,16);
ss_T3TBGT = y0(end,17);
ss_T3TTRT = y0(end,18);
ss_T3ALBT = y0(end,19);

% RB venous blood (last segment of RB blood) 
ss_freeT4RB = y0(end,N+M-1);
ss_T4TBGRB  = y0(end,N+2*M-1);
ss_T4TTRRB  = y0(end,N+3*M-1);
ss_T4ALBRB  = y0(end,N+4*M-1);

ss_freeT3RB = y0(end,N+5*M-1);
ss_T3TBGRB  = y0(end,N+6*M-1);
ss_T3TTRRB  = y0(end,N+7*M-1);
ss_T3ALBRB  = y0(end,N+8*M-1);

% Liver venous blood (last segment of Liver blood) 
ss_freeT4L = y0(end,N+17*M-1);
ss_T4TBGL = y0(end,N+18*M-1);
ss_T4TTRL = y0(end,N+19*M-1);
ss_T4ALBL = y0(end,N+20*M-1);

ss_freeT3L = y0(end,N+21*M-1);
ss_T3TBGL = y0(end,N+22*M-1);
ss_T3TTRL = y0(end,N+23*M-1);
ss_T3ALBL = y0(end,N+24*M-1);

% Percent Change
ss_CACVfT4T = (ss_freeT4T - ss_freeT4B); % /ss_freeT4B * 100;
ss_CACVfT4L = (ss_freeT4L - ss_freeT4B); %/ss_freeT4B * 100;
ss_CACVfT4RB = (ss_freeT4RB - ss_freeT4B); %/ss_freeT4B * 100;

ss_CACVfT3T = (ss_freeT3T - ss_freeT3B); %/ss_freeT3B * 100;
ss_CACVfT3L = (ss_freeT3L - ss_freeT3B); %/ss_freeT3B * 100;
ss_CACVfT3RB = (ss_freeT3RB - ss_freeT3B); %/ss_freeT3B * 100;

ss_CACVT4TBGT = (ss_T4TBGT - ss_T4TBGB); %/ss_T4TBGB * 100;
ss_CACVT4TBGL = (ss_T4TBGL - ss_T4TBGB); %/ss_T4TBGB * 100;
ss_CACVT4TBGRB = (ss_T4TBGRB - ss_T4TBGB); %/ss_T4TBGB * 100;

ss_CACVT3TBGT = (ss_T3TBGT - ss_T3TBGB); %/ss_T3TBGB * 100;
ss_CACVT3TBGL = (ss_T3TBGL - ss_T3TBGB); %/ss_T3TBGB * 100;
ss_CACVT3TBGRB = (ss_T3TBGRB - ss_T3TBGB); %/ss_T3TBGB * 100;

ss_CACVT4TTRT = (ss_T4TTRT - ss_T4TTRB); %/ss_T4TTRB * 100;
ss_CACVT4TTRL = (ss_T4TTRL - ss_T4TTRB); %/ss_T4TTRB * 100;
ss_CACVT4TTRRB = (ss_T4TTRRB - ss_T4TTRB); %/ss_T4TTRB * 100;

ss_CACVT3TTRT = (ss_T3TTRT - ss_T3TTRB); %/ss_T3TTRB * 100;
ss_CACVT3TTRL = (ss_T3TTRL - ss_T3TTRB); %/ss_T3TTRB * 100;
ss_CACVT3TTRRB = (ss_T3TTRRB - ss_T3TTRB); %/ss_T3TTRB * 100;

ss_CACVT4ALBT = (ss_T4ALBT - ss_T4ALBB); %/ss_T4ALBB * 100;
ss_CACVT4ALBL = (ss_T4ALBL - ss_T4ALBB); %/ss_T4ALBB * 100;
ss_CACVT4ALBRB = (ss_T4ALBRB - ss_T4ALBB); %/ss_T4ALBB * 100;

ss_CACVT3ALBT = (ss_T3ALBT - ss_T3ALBB); %/ss_T3ALBB * 100;
ss_CACVT3ALBL = (ss_T3ALBL - ss_T3ALBB); %/ss_T3ALBB * 100;
ss_CACVT3ALBRB = (ss_T3ALBRB - ss_T3ALBB); %/ss_T3ALBB * 100;

ss_CACVtotalT4T = ss_CACVfT4T + ss_CACVT4TBGT + ss_CACVT4TTRT + ss_CACVT4ALBT;
ss_contribution_ratio_T4TBGT = ss_CACVT4TBGT / ss_CACVtotalT4T;
ss_contribution_ratio_T4TTRT = ss_CACVT4TTRT / ss_CACVtotalT4T;
ss_contribution_ratio_T4ALBT = ss_CACVT4ALBT / ss_CACVtotalT4T;

ss_CACVtotalT4L = ss_CACVfT4L + ss_CACVT4TBGL + ss_CACVT4TTRL + ss_CACVT4ALBL;
ss_contribution_ratio_T4TBGL = ss_CACVT4TBGL / ss_CACVtotalT4L;
ss_contribution_ratio_T4TTRL = ss_CACVT4TTRL / ss_CACVtotalT4L;
ss_contribution_ratio_T4ALBL = ss_CACVT4ALBL / ss_CACVtotalT4L;

ss_CACVtotalT4RB = ss_CACVfT4RB + ss_CACVT4TBGRB + ss_CACVT4TTRRB + ss_CACVT4ALBRB;
ss_contribution_ratio_T4TBGRB = ss_CACVT4TBGRB / ss_CACVtotalT4RB;
ss_contribution_ratio_T4TTRRB = ss_CACVT4TTRRB / ss_CACVtotalT4RB;
ss_contribution_ratio_T4ALBRB = ss_CACVT4ALBRB / ss_CACVtotalT4RB;

ss_CACVtotalT3T = ss_CACVfT3T + ss_CACVT3TBGT + ss_CACVT3TTRT + ss_CACVT3ALBT;
ss_contribution_ratio_T3TBGT = ss_CACVT3TBGT / ss_CACVtotalT3T;
ss_contribution_ratio_T3TTRT = ss_CACVT3TTRT / ss_CACVtotalT3T;
ss_contribution_ratio_T3ALBT = ss_CACVT3ALBT / ss_CACVtotalT3T;

ss_CACVtotalT3L = ss_CACVfT3L + ss_CACVT3TBGL + ss_CACVT3TTRL + ss_CACVT3ALBL;
ss_contribution_ratio_T3TBGL = ss_CACVT3TBGL / ss_CACVtotalT3L;
ss_contribution_ratio_T3TTRL = ss_CACVT3TTRL / ss_CACVtotalT3L;
ss_contribution_ratio_T3ALBL = ss_CACVT3ALBL / ss_CACVtotalT3L;

ss_CACVtotalT3RB = ss_CACVfT3RB + ss_CACVT3TBGRB + ss_CACVT3TTRRB + ss_CACVT3ALBRB;
ss_contribution_ratio_T3TBGRB = ss_CACVT3TBGRB / ss_CACVtotalT3RB;
ss_contribution_ratio_T3TTRRB = ss_CACVT3TTRRB / ss_CACVtotalT3RB;
ss_contribution_ratio_T3ALBRB = ss_CACVT3ALBRB / ss_CACVtotalT3RB;

for i = 1:length(sensitivity_analysis_params)

    sensitivity_analysis_params(i)

    param = default_param;
    this_param_results = struct();
    
    tspan1 = [0:10000:1000*24*3600]; %Running for 1000 days
    
    percent_change = 0.05;

 %-----------------Increase the parameter value-----------------%
    if sensitivity_analysis_params(i) == "TBGtot" || sensitivity_analysis_params(i) == "TTRtot" || sensitivity_analysis_params(i) == "ALBtot"

        param.(sensitivity_analysis_params(i)) = default_param.(sensitivity_analysis_params(i)) * (1 + percent_change);

        init = init_default;
        init.TBGB    = param.TBGtot/param.VB;
        init.TTRB    = param.TTRtot/param.VB;
        init.ALBB    = param.ALBtot/param.VB;

        y_init = cell2mat(struct2cell(init));

    elseif sensitivity_analysis_params(i) == "QT" || sensitivity_analysis_params(i) == "QL" || sensitivity_analysis_params(i) == "QRB"
        
        param.(sensitivity_analysis_params(i)) = default_param.(sensitivity_analysis_params(i)) * (1 + percent_change);
        param.QC = param.QT + param.QL + param.QRB;
        y_init = y0;

    else

        param.(sensitivity_analysis_params(i)) = default_param.(sensitivity_analysis_params(i)) * (1 + percent_change);
        y_init = y0;

    end
    
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t1,y1] = ode15s('TH_PBK_Spatial_ODE',tspan1, y_init, options, param);
    
    %Obtain steady-state values after parameter change
    % T4 and T3 differnetial between arterial and venous blood concentrations
    this_param_results.positive_ss_freeT4B = y1(end,1);
    this_param_results.positive_ss_T4TBGB = y1(end,2);
    this_param_results.positive_ss_T4TTRB = y1(end,3);
    this_param_results.positive_ss_T4ALBB = y1(end,4);
    
    this_param_results.positive_ss_freeT3B = y1(end,5);
    this_param_results.positive_ss_T3TBGB = y1(end,6);
    this_param_results.positive_ss_T3TTRB = y1(end,7);
    this_param_results.positive_ss_T3ALBB = y1(end,8);
    
    % Thyroid blood Compartment
    this_param_results.positive_ss_freeT4T = y1(end,12);
    this_param_results.positive_ss_T4TBGT = y1(end,13);
    this_param_results.positive_ss_T4TTRT = y1(end,14);
    this_param_results.positive_ss_T4ALBT = y1(end,15);
    
    this_param_results.positive_ss_freeT3T = y1(end,16);
    this_param_results.positive_ss_T3TBGT = y1(end,17);
    this_param_results.positive_ss_T3TTRT = y1(end,18);
    this_param_results.positive_ss_T3ALBT = y1(end,19);
    
    % RB venous blood (last segment of RB blood) 
    this_param_results.positive_ss_freeT4RB = y1(end,N+M-1);
    this_param_results.positive_ss_T4TBGRB  = y1(end,N+2*M-1);
    this_param_results.positive_ss_T4TTRRB  = y1(end,N+3*M-1);
    this_param_results.positive_ss_T4ALBRB  = y1(end,N+4*M-1);
    
    this_param_results.positive_ss_freeT3RB = y1(end,N+5*M-1);
    this_param_results.positive_ss_T3TBGRB  = y1(end,N+6*M-1);
    this_param_results.positive_ss_T3TTRRB  = y1(end,N+7*M-1);
    this_param_results.positive_ss_T3ALBRB  = y1(end,N+8*M-1);
    
    % Liver venous blood (last segment of Liver blood) 
    this_param_results.positive_ss_freeT4L = y1(end,N+17*M-1);
    this_param_results.positive_ss_T4TBGL = y1(end,N+18*M-1);
    this_param_results.positive_ss_T4TTRL = y1(end,N+19*M-1);
    this_param_results.positive_ss_T4ALBL = y1(end,N+20*M-1);
    
    this_param_results.positive_ss_freeT3L = y1(end,N+21*M-1);
    this_param_results.positive_ss_T3TBGL = y1(end,N+22*M-1);
    this_param_results.positive_ss_T3TTRL = y1(end,N+23*M-1);
    this_param_results.positive_ss_T3ALBL = y1(end,N+24*M-1);
    
    % Percent Change
    this_param_results.positive_ss_CACVfT4T = (this_param_results.positive_ss_freeT4T - this_param_results.positive_ss_freeT4B); %/this_param_results.positive_ss_freeT4B * 100;
    this_param_results.positive_ss_CACVfT4L = (this_param_results.positive_ss_freeT4L - this_param_results.positive_ss_freeT4B); %/this_param_results.positive_ss_freeT4B * 100;
    this_param_results.positive_ss_CACVfT4RB = (this_param_results.positive_ss_freeT4RB - this_param_results.positive_ss_freeT4B); %/this_param_results.positive_ss_freeT4B * 100;
    
    this_param_results.positive_ss_CACVfT3T = (this_param_results.positive_ss_freeT3T - this_param_results.positive_ss_freeT3B); %/this_param_results.positive_ss_freeT3B * 100;
    this_param_results.positive_ss_CACVfT3L = (this_param_results.positive_ss_freeT3L - this_param_results.positive_ss_freeT3B); %/this_param_results.positive_ss_freeT3B * 100;
    this_param_results.positive_ss_CACVfT3RB = (this_param_results.positive_ss_freeT3RB - this_param_results.positive_ss_freeT3B); %/this_param_results.positive_ss_freeT3B * 100;
    
    this_param_results.positive_ss_CACVT4TBGT = (this_param_results.positive_ss_T4TBGT - this_param_results.positive_ss_T4TBGB); %/this_param_results.positive_ss_T4TBGB * 100;
    this_param_results.positive_ss_CACVT4TBGL = (this_param_results.positive_ss_T4TBGL - this_param_results.positive_ss_T4TBGB); %/this_param_results.positive_ss_T4TBGB * 100;
    this_param_results.positive_ss_CACVT4TBGRB = (this_param_results.positive_ss_T4TBGRB - this_param_results.positive_ss_T4TBGB); %/this_param_results.positive_ss_T4TBGB * 100;
    
    this_param_results.positive_ss_CACVT3TBGT = (this_param_results.positive_ss_T3TBGT - this_param_results.positive_ss_T3TBGB); %/this_param_results.positive_ss_T3TBGB * 100;
    this_param_results.positive_ss_CACVT3TBGL = (this_param_results.positive_ss_T3TBGL - this_param_results.positive_ss_T3TBGB); %/this_param_results.positive_ss_T3TBGB * 100;
    this_param_results.positive_ss_CACVT3TBGRB = (this_param_results.positive_ss_T3TBGRB - this_param_results.positive_ss_T3TBGB); %/this_param_results.positive_ss_T3TBGB * 100;
    
    this_param_results.positive_ss_CACVT4TTRT = (this_param_results.positive_ss_T4TTRT - this_param_results.positive_ss_T4TTRB); %/this_param_results.positive_ss_T4TTRB * 100;
    this_param_results.positive_ss_CACVT4TTRL = (this_param_results.positive_ss_T4TTRL - this_param_results.positive_ss_T4TTRB); %/this_param_results.positive_ss_T4TTRB * 100;
    this_param_results.positive_ss_CACVT4TTRRB = (this_param_results.positive_ss_T4TTRRB - this_param_results.positive_ss_T4TTRB); %/this_param_results.positive_ss_T4TTRB * 100;
    
    this_param_results.positive_ss_CACVT3TTRT = (this_param_results.positive_ss_T3TTRT - this_param_results.positive_ss_T3TTRB); %/this_param_results.positive_ss_T3TTRB * 100;
    this_param_results.positive_ss_CACVT3TTRL = (this_param_results.positive_ss_T3TTRL - this_param_results.positive_ss_T3TTRB); %/this_param_results.positive_ss_T3TTRB * 100;
    this_param_results.positive_ss_CACVT3TTRRB = (this_param_results.positive_ss_T3TTRRB - this_param_results.positive_ss_T3TTRB); %/this_param_results.positive_ss_T3TTRB * 100;
    
    this_param_results.positive_ss_CACVT4ALBT = (this_param_results.positive_ss_T4ALBT - this_param_results.positive_ss_T4ALBB); %/this_param_results.positive_ss_T4ALBB * 100;
    this_param_results.positive_ss_CACVT4ALBL = (this_param_results.positive_ss_T4ALBL - this_param_results.positive_ss_T4ALBB); %/this_param_results.positive_ss_T4ALBB * 100;
    this_param_results.positive_ss_CACVT4ALBRB = (this_param_results.positive_ss_T4ALBRB - this_param_results.positive_ss_T4ALBB); %/this_param_results.positive_ss_T4ALBB * 100;
    
    this_param_results.positive_ss_CACVT3ALBT = (this_param_results.positive_ss_T3ALBT - this_param_results.positive_ss_T3ALBB); %/this_param_results.positive_ss_T3ALBB * 100;
    this_param_results.positive_ss_CACVT3ALBL = (this_param_results.positive_ss_T3ALBL - this_param_results.positive_ss_T3ALBB); %/this_param_results.positive_ss_T3ALBB * 100;
    this_param_results.positive_ss_CACVT3ALBRB = (this_param_results.positive_ss_T3ALBRB - this_param_results.positive_ss_T3ALBB); %/this_param_results.positive_ss_T3ALBB * 100;
    
    % Total and relative contributions
    this_param_results.positive_ss_CACVtotalT4T = this_param_results.positive_ss_CACVfT4T + this_param_results.positive_ss_CACVT4TBGT + this_param_results.positive_ss_CACVT4TTRT + this_param_results.positive_ss_CACVT4ALBT;
    this_param_results.positive_ss_contribution_ratio_T4TBGT = this_param_results.positive_ss_CACVT4TBGT / this_param_results.positive_ss_CACVtotalT4T;
    this_param_results.positive_ss_contribution_ratio_T4TTRT = this_param_results.positive_ss_CACVT4TTRT / this_param_results.positive_ss_CACVtotalT4T;
    this_param_results.positive_ss_contribution_ratio_T4ALBT = this_param_results.positive_ss_CACVT4ALBT / this_param_results.positive_ss_CACVtotalT4T;

    this_param_results.positive_ss_CACVtotalT4L = this_param_results.positive_ss_CACVfT4L + this_param_results.positive_ss_CACVT4TBGL + this_param_results.positive_ss_CACVT4TTRL + this_param_results.positive_ss_CACVT4ALBL;
    this_param_results.positive_ss_contribution_ratio_T4TBGL = this_param_results.positive_ss_CACVT4TBGL / this_param_results.positive_ss_CACVtotalT4L;
    this_param_results.positive_ss_contribution_ratio_T4TTRL = this_param_results.positive_ss_CACVT4TTRL / this_param_results.positive_ss_CACVtotalT4L;
    this_param_results.positive_ss_contribution_ratio_T4ALBL = this_param_results.positive_ss_CACVT4ALBL / this_param_results.positive_ss_CACVtotalT4L;

    this_param_results.positive_ss_CACVtotalT4RB = this_param_results.positive_ss_CACVfT4RB + this_param_results.positive_ss_CACVT4TBGRB + this_param_results.positive_ss_CACVT4TTRRB + this_param_results.positive_ss_CACVT4ALBRB;
    this_param_results.positive_ss_contribution_ratio_T4TBGRB = this_param_results.positive_ss_CACVT4TBGRB / this_param_results.positive_ss_CACVtotalT4RB;
    this_param_results.positive_ss_contribution_ratio_T4TTRRB = this_param_results.positive_ss_CACVT4TTRRB / this_param_results.positive_ss_CACVtotalT4RB;
    this_param_results.positive_ss_contribution_ratio_T4ALBRB = this_param_results.positive_ss_CACVT4ALBRB / this_param_results.positive_ss_CACVtotalT4RB;
    
    this_param_results.positive_ss_CACVtotalT3T = this_param_results.positive_ss_CACVfT3T + this_param_results.positive_ss_CACVT3TBGT + this_param_results.positive_ss_CACVT3TTRT + this_param_results.positive_ss_CACVT3ALBT;
    this_param_results.positive_ss_contribution_ratio_T3TBGT = this_param_results.positive_ss_CACVT3TBGT / this_param_results.positive_ss_CACVtotalT3T;
    this_param_results.positive_ss_contribution_ratio_T3TTRT = this_param_results.positive_ss_CACVT3TTRT / this_param_results.positive_ss_CACVtotalT3T;
    this_param_results.positive_ss_contribution_ratio_T3ALBT = this_param_results.positive_ss_CACVT3ALBT / this_param_results.positive_ss_CACVtotalT3T;

    this_param_results.positive_ss_CACVtotalT3L = this_param_results.positive_ss_CACVfT3L + this_param_results.positive_ss_CACVT3TBGL + this_param_results.positive_ss_CACVT3TTRL + this_param_results.positive_ss_CACVT3ALBL;
    this_param_results.positive_ss_contribution_ratio_T3TBGL = this_param_results.positive_ss_CACVT3TBGL / this_param_results.positive_ss_CACVtotalT3L;
    this_param_results.positive_ss_contribution_ratio_T3TTRL = this_param_results.positive_ss_CACVT3TTRL / this_param_results.positive_ss_CACVtotalT3L;
    this_param_results.positive_ss_contribution_ratio_T3ALBL = this_param_results.positive_ss_CACVT3ALBL / this_param_results.positive_ss_CACVtotalT3L;

    this_param_results.positive_ss_CACVtotalT3RB = this_param_results.positive_ss_CACVfT3RB + this_param_results.positive_ss_CACVT3TBGRB + this_param_results.positive_ss_CACVT3TTRRB + this_param_results.positive_ss_CACVT3ALBRB;
    this_param_results.positive_ss_contribution_ratio_T3TBGRB = this_param_results.positive_ss_CACVT3TBGRB / this_param_results.positive_ss_CACVtotalT3RB;
    this_param_results.positive_ss_contribution_ratio_T3TTRRB = this_param_results.positive_ss_CACVT3TTRRB / this_param_results.positive_ss_CACVtotalT3RB;
    this_param_results.positive_ss_contribution_ratio_T3ALBRB = this_param_results.positive_ss_CACVT3ALBRB / this_param_results.positive_ss_CACVtotalT3RB;




 %-----------------Decrease the parameter value-----------------%

    if sensitivity_analysis_params(i) == "TBGtot" || sensitivity_analysis_params(i) == "TTRtot" || sensitivity_analysis_params(i) == "ALBtot"

        param.(sensitivity_analysis_params(i)) = default_param.(sensitivity_analysis_params(i)) * (1 - percent_change);

        init = init_default;
        init.TBGB    = param.TBGtot/param.VB;
        init.TTRB    = param.TTRtot/param.VB;
        init.ALBB    = param.ALBtot/param.VB;

        y_init = cell2mat(struct2cell(init));

    elseif sensitivity_analysis_params(i) == "QT" || sensitivity_analysis_params(i) == "QL" || sensitivity_analysis_params(i) == "QRB"
        
        param.(sensitivity_analysis_params(i)) = default_param.(sensitivity_analysis_params(i)) * (1 - percent_change);
        param.QC = param.QT + param.QL + param.QRB;
        y_init = y0;

    else

        param.(sensitivity_analysis_params(i)) = default_param.(sensitivity_analysis_params(i)) * (1 - percent_change);
        y_init = y0;

    end

    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t1,y1] = ode15s('TH_PBK_Spatial_ODE',tspan1, y_init, options, param);

    %Obtain steady-state values after parameter change
    % T4 and T3 differnetial between arterial and venous blood concentrations
    this_param_results.negative_ss_freeT4B = y1(end,1);
    this_param_results.negative_ss_T4TBGB = y1(end,2);
    this_param_results.negative_ss_T4TTRB = y1(end,3);
    this_param_results.negative_ss_T4ALBB = y1(end,4);
    
    this_param_results.negative_ss_freeT3B = y1(end,5);
    this_param_results.negative_ss_T3TBGB = y1(end,6);
    this_param_results.negative_ss_T3TTRB = y1(end,7);
    this_param_results.negative_ss_T3ALBB = y1(end,8);
    
    % Thyroid blood Compartment
    this_param_results.negative_ss_freeT4T = y1(end,12);
    this_param_results.negative_ss_T4TBGT = y1(end,13);
    this_param_results.negative_ss_T4TTRT = y1(end,14);
    this_param_results.negative_ss_T4ALBT = y1(end,15);
    
    this_param_results.negative_ss_freeT3T = y1(end,16);
    this_param_results.negative_ss_T3TBGT = y1(end,17);
    this_param_results.negative_ss_T3TTRT = y1(end,18);
    this_param_results.negative_ss_T3ALBT = y1(end,19);
    
    % RB venous blood (last segment of RB blood) 
    this_param_results.negative_ss_freeT4RB = y1(end,N+M-1);
    this_param_results.negative_ss_T4TBGRB  = y1(end,N+2*M-1);
    this_param_results.negative_ss_T4TTRRB  = y1(end,N+3*M-1);
    this_param_results.negative_ss_T4ALBRB  = y1(end,N+4*M-1);
    
    this_param_results.negative_ss_freeT3RB = y1(end,N+5*M-1);
    this_param_results.negative_ss_T3TBGRB  = y1(end,N+6*M-1);
    this_param_results.negative_ss_T3TTRRB  = y1(end,N+7*M-1);
    this_param_results.negative_ss_T3ALBRB  = y1(end,N+8*M-1);
    
    % Liver venous blood (last segment of Liver blood) 
    this_param_results.negative_ss_freeT4L = y1(end,N+17*M-1);
    this_param_results.negative_ss_T4TBGL = y1(end,N+18*M-1);
    this_param_results.negative_ss_T4TTRL = y1(end,N+19*M-1);
    this_param_results.negative_ss_T4ALBL = y1(end,N+20*M-1);
    
    this_param_results.negative_ss_freeT3L = y1(end,N+21*M-1);
    this_param_results.negative_ss_T3TBGL = y1(end,N+22*M-1);
    this_param_results.negative_ss_T3TTRL = y1(end,N+23*M-1);
    this_param_results.negative_ss_T3ALBL = y1(end,N+24*M-1);
    
    % Percent Change
    this_param_results.negative_ss_CACVfT4T = (this_param_results.negative_ss_freeT4T - this_param_results.negative_ss_freeT4B); %/this_param_results.negative_ss_freeT4B * 100;
    this_param_results.negative_ss_CACVfT4L = (this_param_results.negative_ss_freeT4L - this_param_results.negative_ss_freeT4B); %/this_param_results.negative_ss_freeT4B * 100;
    this_param_results.negative_ss_CACVfT4RB = (this_param_results.negative_ss_freeT4RB - this_param_results.negative_ss_freeT4B); %/this_param_results.negative_ss_freeT4B * 100;
    
    this_param_results.negative_ss_CACVfT3T = (this_param_results.negative_ss_freeT3T - this_param_results.negative_ss_freeT3B); %/this_param_results.negative_ss_freeT3B * 100;
    this_param_results.negative_ss_CACVfT3L = (this_param_results.negative_ss_freeT3L - this_param_results.negative_ss_freeT3B); %/this_param_results.negative_ss_freeT3B * 100;
    this_param_results.negative_ss_CACVfT3RB = (this_param_results.negative_ss_freeT3RB - this_param_results.negative_ss_freeT3B); %/this_param_results.negative_ss_freeT3B * 100;
    
    this_param_results.negative_ss_CACVT4TBGT = (this_param_results.negative_ss_T4TBGT - this_param_results.negative_ss_T4TBGB); %/this_param_results.negative_ss_T4TBGB * 100;
    this_param_results.negative_ss_CACVT4TBGL = (this_param_results.negative_ss_T4TBGL - this_param_results.negative_ss_T4TBGB); %/this_param_results.negative_ss_T4TBGB * 100;
    this_param_results.negative_ss_CACVT4TBGRB = (this_param_results.negative_ss_T4TBGRB - this_param_results.negative_ss_T4TBGB); %/this_param_results.negative_ss_T4TBGB * 100;
    
    this_param_results.negative_ss_CACVT3TBGT = (this_param_results.negative_ss_T3TBGT - this_param_results.negative_ss_T3TBGB); %/this_param_results.negative_ss_T3TBGB * 100;
    this_param_results.negative_ss_CACVT3TBGL = (this_param_results.negative_ss_T3TBGL - this_param_results.negative_ss_T3TBGB); %/this_param_results.negative_ss_T3TBGB * 100;
    this_param_results.negative_ss_CACVT3TBGRB = (this_param_results.negative_ss_T3TBGRB - this_param_results.negative_ss_T3TBGB); %/this_param_results.negative_ss_T3TBGB * 100;
    
    this_param_results.negative_ss_CACVT4TTRT = (this_param_results.negative_ss_T4TTRT - this_param_results.negative_ss_T4TTRB); %/this_param_results.negative_ss_T4TTRB * 100;
    this_param_results.negative_ss_CACVT4TTRL = (this_param_results.negative_ss_T4TTRL - this_param_results.negative_ss_T4TTRB); %/this_param_results.negative_ss_T4TTRB * 100;
    this_param_results.negative_ss_CACVT4TTRRB = (this_param_results.negative_ss_T4TTRRB - this_param_results.negative_ss_T4TTRB); %/this_param_results.negative_ss_T4TTRB * 100;
    
    this_param_results.negative_ss_CACVT3TTRT = (this_param_results.negative_ss_T3TTRT - this_param_results.negative_ss_T3TTRB); %/this_param_results.negative_ss_T3TTRB * 100;
    this_param_results.negative_ss_CACVT3TTRL = (this_param_results.negative_ss_T3TTRL - this_param_results.negative_ss_T3TTRB); %/this_param_results.negative_ss_T3TTRB * 100;
    this_param_results.negative_ss_CACVT3TTRRB = (this_param_results.negative_ss_T3TTRRB - this_param_results.negative_ss_T3TTRB); %/this_param_results.negative_ss_T3TTRB * 100;
    
    this_param_results.negative_ss_CACVT4ALBT = (this_param_results.negative_ss_T4ALBT - this_param_results.negative_ss_T4ALBB); %; %/this_param_results.negative_ss_T4ALBB * 100;
    this_param_results.negative_ss_CACVT4ALBL = (this_param_results.negative_ss_T4ALBL - this_param_results.negative_ss_T4ALBB); %/this_param_results.negative_ss_T4ALBB * 100;
    this_param_results.negative_ss_CACVT4ALBRB = (this_param_results.negative_ss_T4ALBRB - this_param_results.negative_ss_T4ALBB); %/this_param_results.negative_ss_T4ALBB * 100;
    
    this_param_results.negative_ss_CACVT3ALBT = (this_param_results.negative_ss_T3ALBT - this_param_results.negative_ss_T3ALBB); %/this_param_results.negative_ss_T3ALBB * 100;
    this_param_results.negative_ss_CACVT3ALBL = (this_param_results.negative_ss_T3ALBL - this_param_results.negative_ss_T3ALBB); %/this_param_results.negative_ss_T3ALBB * 100;
    this_param_results.negative_ss_CACVT3ALBRB = (this_param_results.negative_ss_T3ALBRB - this_param_results.negative_ss_T3ALBB); %/this_param_results.negative_ss_T3ALBB * 100;
    
    this_param_results.negative_ss_CACVtotalT4T = this_param_results.negative_ss_CACVfT4T + this_param_results.negative_ss_CACVT4TBGT + this_param_results.negative_ss_CACVT4TTRT + this_param_results.negative_ss_CACVT4ALBT;
    this_param_results.negative_ss_contribution_ratio_T4TBGT = this_param_results.negative_ss_CACVT4TBGT / this_param_results.negative_ss_CACVtotalT4T;
    this_param_results.negative_ss_contribution_ratio_T4TTRT = this_param_results.negative_ss_CACVT4TTRT / this_param_results.negative_ss_CACVtotalT4T;
    this_param_results.negative_ss_contribution_ratio_T4ALBT = this_param_results.negative_ss_CACVT4ALBT / this_param_results.negative_ss_CACVtotalT4T;
    
    this_param_results.negative_ss_CACVtotalT4L = this_param_results.negative_ss_CACVfT4L + this_param_results.negative_ss_CACVT4TBGL + this_param_results.negative_ss_CACVT4TTRL + this_param_results.negative_ss_CACVT4ALBL;
    this_param_results.negative_ss_contribution_ratio_T4TBGL = this_param_results.negative_ss_CACVT4TBGL / this_param_results.negative_ss_CACVtotalT4L;
    this_param_results.negative_ss_contribution_ratio_T4TTRL = this_param_results.negative_ss_CACVT4TTRL / this_param_results.negative_ss_CACVtotalT4L;
    this_param_results.negative_ss_contribution_ratio_T4ALBL = this_param_results.negative_ss_CACVT4ALBL / this_param_results.negative_ss_CACVtotalT4L;

    this_param_results.negative_ss_CACVtotalT4RB = this_param_results.negative_ss_CACVfT4RB + this_param_results.negative_ss_CACVT4TBGRB + this_param_results.negative_ss_CACVT4TTRRB + this_param_results.negative_ss_CACVT4ALBRB;
    this_param_results.negative_ss_contribution_ratio_T4TBGRB = this_param_results.negative_ss_CACVT4TBGRB / this_param_results.negative_ss_CACVtotalT4RB;
    this_param_results.negative_ss_contribution_ratio_T4TTRRB = this_param_results.negative_ss_CACVT4TTRRB / this_param_results.negative_ss_CACVtotalT4RB;
    this_param_results.negative_ss_contribution_ratio_T4ALBRB = this_param_results.negative_ss_CACVT4ALBRB / this_param_results.negative_ss_CACVtotalT4RB;

    this_param_results.negative_ss_CACVtotalT3T = this_param_results.negative_ss_CACVfT3T + this_param_results.negative_ss_CACVT3TBGT + this_param_results.negative_ss_CACVT3TTRT + this_param_results.negative_ss_CACVT3ALBT;
    this_param_results.negative_ss_contribution_ratio_T3TBGT = this_param_results.negative_ss_CACVT3TBGT / this_param_results.negative_ss_CACVtotalT3T;
    this_param_results.negative_ss_contribution_ratio_T3TTRT = this_param_results.negative_ss_CACVT3TTRT / this_param_results.negative_ss_CACVtotalT3T;
    this_param_results.negative_ss_contribution_ratio_T3ALBT = this_param_results.negative_ss_CACVT3ALBT / this_param_results.negative_ss_CACVtotalT3T;
    
    this_param_results.negative_ss_CACVtotalT3L = this_param_results.negative_ss_CACVfT3L + this_param_results.negative_ss_CACVT3TBGL + this_param_results.negative_ss_CACVT3TTRL + this_param_results.negative_ss_CACVT3ALBL;
    this_param_results.negative_ss_contribution_ratio_T3TBGL = this_param_results.negative_ss_CACVT3TBGL / this_param_results.negative_ss_CACVtotalT3L;
    this_param_results.negative_ss_contribution_ratio_T3TTRL = this_param_results.negative_ss_CACVT3TTRL / this_param_results.negative_ss_CACVtotalT3L;
    this_param_results.negative_ss_contribution_ratio_T3ALBL = this_param_results.negative_ss_CACVT3ALBL / this_param_results.negative_ss_CACVtotalT3L;

    this_param_results.negative_ss_CACVtotalT3RB = this_param_results.negative_ss_CACVfT3RB + this_param_results.negative_ss_CACVT3TBGRB + this_param_results.negative_ss_CACVT3TTRRB + this_param_results.negative_ss_CACVT3ALBRB;
    this_param_results.negative_ss_contribution_ratio_T3TBGRB = this_param_results.negative_ss_CACVT3TBGRB / this_param_results.negative_ss_CACVtotalT3RB;
    this_param_results.negative_ss_contribution_ratio_T3TTRRB = this_param_results.negative_ss_CACVT3TTRRB / this_param_results.negative_ss_CACVtotalT3RB;
    this_param_results.negative_ss_contribution_ratio_T3ALBRB = this_param_results.negative_ss_CACVT3ALBRB / this_param_results.negative_ss_CACVtotalT3RB;


 %-----------------------Calculate sensitivity coefficient----------------------%
    this_param_results.sensitivity_coeff_CACVfT4T = mean([(this_param_results.positive_ss_CACVfT4T - ss_CACVfT4T)/ss_CACVfT4T/percent_change, (ss_CACVfT4T - this_param_results.negative_ss_CACVfT4T)/ss_CACVfT4T/percent_change]);
    this_param_results.sensitivity_coeff_CACVfT4L = mean([(this_param_results.positive_ss_CACVfT4L - ss_CACVfT4L)/ss_CACVfT4L/percent_change, (ss_CACVfT4L - this_param_results.negative_ss_CACVfT4L)/ss_CACVfT4L/percent_change]);
    this_param_results.sensitivity_coeff_CACVfT4RB = mean([(this_param_results.positive_ss_CACVfT4RB - ss_CACVfT4RB)/ss_CACVfT4RB/percent_change, (ss_CACVfT4RB - this_param_results.negative_ss_CACVfT4RB)/ss_CACVfT4RB/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVfT3T = mean([(this_param_results.positive_ss_CACVfT3T - ss_CACVfT3T)/ss_CACVfT3T/percent_change, (ss_CACVfT3T - this_param_results.negative_ss_CACVfT3T)/ss_CACVfT3T/percent_change]);
    this_param_results.sensitivity_coeff_CACVfT3L = mean([(this_param_results.positive_ss_CACVfT3L - ss_CACVfT3L)/ss_CACVfT3L/percent_change, (ss_CACVfT3L - this_param_results.negative_ss_CACVfT3L)/ss_CACVfT3L/percent_change]);
    this_param_results.sensitivity_coeff_CACVfT3RB = mean([(this_param_results.positive_ss_CACVfT3RB - ss_CACVfT3RB)/ss_CACVfT3RB/percent_change, (ss_CACVfT3RB - this_param_results.negative_ss_CACVfT3RB)/ss_CACVfT3RB/percent_change]);

    this_param_results.sensitivity_coeff_CACVT4TBGT = mean([(this_param_results.positive_ss_CACVT4TBGT - ss_CACVT4TBGT)/ss_CACVT4TBGT/percent_change, (ss_CACVT4TBGT - this_param_results.negative_ss_CACVT4TBGT)/ss_CACVT4TBGT/percent_change]);
    this_param_results.sensitivity_coeff_CACVT4TBGL = mean([(this_param_results.positive_ss_CACVT4TBGL - ss_CACVT4TBGL)/ss_CACVT4TBGL/percent_change, (ss_CACVT4TBGL - this_param_results.negative_ss_CACVT4TBGL)/ss_CACVT4TBGL/percent_change]);
    this_param_results.sensitivity_coeff_CACVT4TBGRB = mean([(this_param_results.positive_ss_CACVT4TBGRB - ss_CACVT4TBGRB)/ss_CACVT4TBGRB/percent_change, (ss_CACVT4TBGRB - this_param_results.negative_ss_CACVT4TBGRB)/ss_CACVT4TBGRB/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVT3TBGT = mean([(this_param_results.positive_ss_CACVT3TBGT - ss_CACVT3TBGT)/ss_CACVT3TBGT/percent_change, (ss_CACVT3TBGT - this_param_results.negative_ss_CACVT3TBGT)/ss_CACVT3TBGT/percent_change]);
    this_param_results.sensitivity_coeff_CACVT3TBGL = mean([(this_param_results.positive_ss_CACVT3TBGL - ss_CACVT3TBGL)/ss_CACVT3TBGL/percent_change, (ss_CACVT3TBGL - this_param_results.negative_ss_CACVT3TBGL)/ss_CACVT3TBGL/percent_change]);
    this_param_results.sensitivity_coeff_CACVT3TBGRB = mean([(this_param_results.positive_ss_CACVT3TBGRB - ss_CACVT3TBGRB)/ss_CACVT3TBGRB/percent_change, (ss_CACVT3TBGRB - this_param_results.negative_ss_CACVT3TBGRB)/ss_CACVT3TBGRB/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVT4TTRT = mean([(this_param_results.positive_ss_CACVT4TTRT - ss_CACVT4TTRT)/ss_CACVT4TTRT/percent_change, (ss_CACVT4TTRT - this_param_results.negative_ss_CACVT4TTRT)/ss_CACVT4TTRT/percent_change]);
    this_param_results.sensitivity_coeff_CACVT4TTRL = mean([(this_param_results.positive_ss_CACVT4TTRL - ss_CACVT4TTRL)/ss_CACVT4TTRL/percent_change, (ss_CACVT4TTRL - this_param_results.negative_ss_CACVT4TTRL)/ss_CACVT4TTRL/percent_change]);
    this_param_results.sensitivity_coeff_CACVT4TTRRB = mean([(this_param_results.positive_ss_CACVT4TTRRB - ss_CACVT4TTRRB)/ss_CACVT4TTRRB/percent_change, (ss_CACVT4TTRRB - this_param_results.negative_ss_CACVT4TTRRB)/ss_CACVT4TTRRB/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVT3TTRT = mean([(this_param_results.positive_ss_CACVT3TTRT - ss_CACVT3TTRT)/ss_CACVT3TTRT/percent_change, (ss_CACVT3TTRT - this_param_results.negative_ss_CACVT3TTRT)/ss_CACVT3TTRT/percent_change]);
    this_param_results.sensitivity_coeff_CACVT3TTRL = mean([(this_param_results.positive_ss_CACVT3TTRL - ss_CACVT3TTRL)/ss_CACVT3TTRL/percent_change, (ss_CACVT3TTRL - this_param_results.negative_ss_CACVT3TTRL)/ss_CACVT3TTRL/percent_change]);
    this_param_results.sensitivity_coeff_CACVT3TTRRB = mean([(this_param_results.positive_ss_CACVT3TTRRB - ss_CACVT3TTRRB)/ss_CACVT3TTRRB/percent_change, (ss_CACVT3TTRRB - this_param_results.negative_ss_CACVT3TTRRB)/ss_CACVT3TTRRB/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVT4ALBT = mean([(this_param_results.positive_ss_CACVT4ALBT - ss_CACVT4ALBT)/ss_CACVT4ALBT/percent_change, (ss_CACVT4ALBT - this_param_results.negative_ss_CACVT4ALBT)/ss_CACVT4ALBT/percent_change]);
    this_param_results.sensitivity_coeff_CACVT4ALBL = mean([(this_param_results.positive_ss_CACVT4ALBL - ss_CACVT4ALBL)/ss_CACVT4ALBL/percent_change, (ss_CACVT4ALBL - this_param_results.negative_ss_CACVT4ALBL)/ss_CACVT4ALBL/percent_change]);
    this_param_results.sensitivity_coeff_CACVT4ALBRB = mean([(this_param_results.positive_ss_CACVT4ALBRB - ss_CACVT4ALBRB)/ss_CACVT4ALBRB/percent_change, (ss_CACVT4ALBRB - this_param_results.negative_ss_CACVT4ALBRB)/ss_CACVT4ALBRB/percent_change]);

    this_param_results.sensitivity_coeff_CACVT3ALBT = mean([(this_param_results.positive_ss_CACVT3ALBT - ss_CACVT3ALBT)/ss_CACVT3ALBT/percent_change, (ss_CACVT3ALBT - this_param_results.negative_ss_CACVT3ALBT)/ss_CACVT3ALBT/percent_change]);
    this_param_results.sensitivity_coeff_CACVT3ALBL = mean([(this_param_results.positive_ss_CACVT3ALBL - ss_CACVT3ALBL)/ss_CACVT3ALBL/percent_change, (ss_CACVT3ALBL - this_param_results.negative_ss_CACVT3ALBL)/ss_CACVT3ALBL/percent_change]);
    this_param_results.sensitivity_coeff_CACVT3ALBRB = mean([(this_param_results.positive_ss_CACVT3ALBRB - ss_CACVT3ALBRB)/ss_CACVT3ALBRB/percent_change, (ss_CACVT3ALBRB - this_param_results.negative_ss_CACVT3ALBRB)/ss_CACVT3ALBRB/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVtotalT4T = mean([(this_param_results.positive_ss_CACVtotalT4T - ss_CACVtotalT4T)/ss_CACVtotalT4T/percent_change, (ss_CACVtotalT4T - this_param_results.negative_ss_CACVtotalT4T)/ss_CACVtotalT4T/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4TBGT = mean([(this_param_results.positive_ss_contribution_ratio_T4TBGT - ss_contribution_ratio_T4TBGT)/ss_contribution_ratio_T4TBGT/percent_change, (ss_contribution_ratio_T4TBGT - this_param_results.negative_ss_contribution_ratio_T4TBGT)/ss_contribution_ratio_T4TBGT/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4TTRT = mean([(this_param_results.positive_ss_contribution_ratio_T4TTRT - ss_contribution_ratio_T4TTRT)/ss_contribution_ratio_T4TTRT/percent_change, (ss_contribution_ratio_T4TTRT - this_param_results.negative_ss_contribution_ratio_T4TTRT)/ss_contribution_ratio_T4TTRT/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4ALBT = mean([(this_param_results.positive_ss_contribution_ratio_T4ALBT - ss_contribution_ratio_T4ALBT)/ss_contribution_ratio_T4ALBT/percent_change, (ss_contribution_ratio_T4ALBT - this_param_results.negative_ss_contribution_ratio_T4ALBT)/ss_contribution_ratio_T4ALBT/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVtotalT4L = mean([(this_param_results.positive_ss_CACVtotalT4L - ss_CACVtotalT4L)/ss_CACVtotalT4L/percent_change, (ss_CACVtotalT4L - this_param_results.negative_ss_CACVtotalT4L)/ss_CACVtotalT4L/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4TBGL = mean([(this_param_results.positive_ss_contribution_ratio_T4TBGL - ss_contribution_ratio_T4TBGL)/ss_contribution_ratio_T4TBGL/percent_change, (ss_contribution_ratio_T4TBGL - this_param_results.negative_ss_contribution_ratio_T4TBGL)/ss_contribution_ratio_T4TBGL/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4TTRL = mean([(this_param_results.positive_ss_contribution_ratio_T4TTRL - ss_contribution_ratio_T4TTRL)/ss_contribution_ratio_T4TTRL/percent_change, (ss_contribution_ratio_T4TTRL - this_param_results.negative_ss_contribution_ratio_T4TTRL)/ss_contribution_ratio_T4TTRL/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4ALBL = mean([(this_param_results.positive_ss_contribution_ratio_T4ALBL - ss_contribution_ratio_T4ALBL)/ss_contribution_ratio_T4ALBL/percent_change, (ss_contribution_ratio_T4ALBL - this_param_results.negative_ss_contribution_ratio_T4ALBL)/ss_contribution_ratio_T4ALBL/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVtotalT4RB = mean([(this_param_results.positive_ss_CACVtotalT4RB - ss_CACVtotalT4RB)/ss_CACVtotalT4RB/percent_change, (ss_CACVtotalT4RB - this_param_results.negative_ss_CACVtotalT4RB)/ss_CACVtotalT4RB/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4TBGRB = mean([(this_param_results.positive_ss_contribution_ratio_T4TBGRB - ss_contribution_ratio_T4TBGRB)/ss_contribution_ratio_T4TBGRB/percent_change, (ss_contribution_ratio_T4TBGRB - this_param_results.negative_ss_contribution_ratio_T4TBGRB)/ss_contribution_ratio_T4TBGRB/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4TTRRB = mean([(this_param_results.positive_ss_contribution_ratio_T4TTRRB - ss_contribution_ratio_T4TTRRB)/ss_contribution_ratio_T4TTRRB/percent_change, (ss_contribution_ratio_T4TTRRB - this_param_results.negative_ss_contribution_ratio_T4TTRRB)/ss_contribution_ratio_T4TTRRB/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T4ALBRB = mean([(this_param_results.positive_ss_contribution_ratio_T4ALBRB - ss_contribution_ratio_T4ALBRB)/ss_contribution_ratio_T4ALBRB/percent_change, (ss_contribution_ratio_T4ALBRB - this_param_results.negative_ss_contribution_ratio_T4ALBRB)/ss_contribution_ratio_T4ALBRB/percent_change]);

    this_param_results.sensitivity_coeff_CACVtotalT3T = mean([(this_param_results.positive_ss_CACVtotalT3T - ss_CACVtotalT3T)/ss_CACVtotalT3T/percent_change, (ss_CACVtotalT3T - this_param_results.negative_ss_CACVtotalT3T)/ss_CACVtotalT3T/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3TBGT = mean([(this_param_results.positive_ss_contribution_ratio_T3TBGT - ss_contribution_ratio_T3TBGT)/ss_contribution_ratio_T3TBGT/percent_change, (ss_contribution_ratio_T3TBGT - this_param_results.negative_ss_contribution_ratio_T3TBGT)/ss_contribution_ratio_T3TBGT/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3TTRT = mean([(this_param_results.positive_ss_contribution_ratio_T3TTRT - ss_contribution_ratio_T3TTRT)/ss_contribution_ratio_T3TTRT/percent_change, (ss_contribution_ratio_T3TTRT - this_param_results.negative_ss_contribution_ratio_T3TTRT)/ss_contribution_ratio_T3TTRT/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3ALBT = mean([(this_param_results.positive_ss_contribution_ratio_T3ALBT - ss_contribution_ratio_T3ALBT)/ss_contribution_ratio_T3ALBT/percent_change, (ss_contribution_ratio_T3ALBT - this_param_results.negative_ss_contribution_ratio_T3ALBT)/ss_contribution_ratio_T3ALBT/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVtotalT3L = mean([(this_param_results.positive_ss_CACVtotalT3L - ss_CACVtotalT3L)/ss_CACVtotalT3L/percent_change, (ss_CACVtotalT3L - this_param_results.negative_ss_CACVtotalT3L)/ss_CACVtotalT3L/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3TBGL = mean([(this_param_results.positive_ss_contribution_ratio_T3TBGL - ss_contribution_ratio_T3TBGL)/ss_contribution_ratio_T3TBGL/percent_change, (ss_contribution_ratio_T3TBGL - this_param_results.negative_ss_contribution_ratio_T3TBGL)/ss_contribution_ratio_T3TBGL/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3TTRL = mean([(this_param_results.positive_ss_contribution_ratio_T3TTRL - ss_contribution_ratio_T3TTRL)/ss_contribution_ratio_T3TTRL/percent_change, (ss_contribution_ratio_T3TTRL - this_param_results.negative_ss_contribution_ratio_T3TTRL)/ss_contribution_ratio_T3TTRL/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3ALBL = mean([(this_param_results.positive_ss_contribution_ratio_T3ALBL - ss_contribution_ratio_T3ALBL)/ss_contribution_ratio_T3ALBL/percent_change, (ss_contribution_ratio_T3ALBL - this_param_results.negative_ss_contribution_ratio_T3ALBL)/ss_contribution_ratio_T3ALBL/percent_change]);
    
    this_param_results.sensitivity_coeff_CACVtotalT3RB = mean([(this_param_results.positive_ss_CACVtotalT3RB - ss_CACVtotalT3RB)/ss_CACVtotalT3RB/percent_change, (ss_CACVtotalT3RB - this_param_results.negative_ss_CACVtotalT3RB)/ss_CACVtotalT3RB/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3TBGRB = mean([(this_param_results.positive_ss_contribution_ratio_T3TBGRB - ss_contribution_ratio_T3TBGRB)/ss_contribution_ratio_T3TBGRB/percent_change, (ss_contribution_ratio_T3TBGRB - this_param_results.negative_ss_contribution_ratio_T3TBGRB)/ss_contribution_ratio_T3TBGRB/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3TTRRB = mean([(this_param_results.positive_ss_contribution_ratio_T3TTRRB - ss_contribution_ratio_T3TTRRB)/ss_contribution_ratio_T3TTRRB/percent_change, (ss_contribution_ratio_T3TTRRB - this_param_results.negative_ss_contribution_ratio_T3TTRRB)/ss_contribution_ratio_T3TTRRB/percent_change]);
    this_param_results.sensitivity_coeff_contribution_ratio_T3ALBRB = mean([(this_param_results.positive_ss_contribution_ratio_T3ALBRB - ss_contribution_ratio_T3ALBRB)/ss_contribution_ratio_T3ALBRB/percent_change, (ss_contribution_ratio_T3ALBRB - this_param_results.negative_ss_contribution_ratio_T3ALBRB)/ss_contribution_ratio_T3ALBRB/percent_change]);

    sensitivity_analysis_results.(sensitivity_analysis_params(i)) = this_param_results;

   % save('sensitivity_analysis_results.mat', 'sensitivity_analysis_results');
end



coeffs = [
    "sensitivity_coeff_contribution_ratio_T4TBGT"
    "sensitivity_coeff_contribution_ratio_T4TTRT"
    "sensitivity_coeff_contribution_ratio_T4ALBT"
    "sensitivity_coeff_contribution_ratio_T4TBGL"
    "sensitivity_coeff_contribution_ratio_T4TTRL"
    "sensitivity_coeff_contribution_ratio_T4ALBL"
    "sensitivity_coeff_contribution_ratio_T4TBGRB"
    "sensitivity_coeff_contribution_ratio_T4TTRRB"
    "sensitivity_coeff_contribution_ratio_T4ALBRB"
    "sensitivity_coeff_contribution_ratio_T3TBGT"
    "sensitivity_coeff_contribution_ratio_T3TTRT"
    "sensitivity_coeff_contribution_ratio_T3ALBT"
    "sensitivity_coeff_contribution_ratio_T3TBGL"
    "sensitivity_coeff_contribution_ratio_T3TTRL"
    "sensitivity_coeff_contribution_ratio_T3ALBL"
    "sensitivity_coeff_contribution_ratio_T3TBGRB"
    "sensitivity_coeff_contribution_ratio_T3TTRRB"
    "sensitivity_coeff_contribution_ratio_T3ALBRB"
];

for coeff = 1:length(coeffs)

    coeff_values = [];

    for param = 1:length(sensitivity_analysis_params)

        param_struct = sensitivity_analysis_results.(sensitivity_analysis_params(param));
        coeff_values = [coeff_values, param_struct.(coeffs(coeff))];

    end

    if contains(coeffs(coeff), 'TBG')
        color = [0.93 0.69 0.13];
    elseif contains(coeffs(coeff), 'TTR')
        color = [0.49,0.18,0.56];
    elseif contains(coeffs(coeff), 'ALB')
        color = [0.00 0.45 0.74];
    else
        color = [0.0 0.0 0.0]; % Default if none of the conditions match
    end

    % Tornado plot;
    [~,idx] = sort(abs(coeff_values), 'ascend'); % Obtain index after sorting by absolute values; 
    param_names = sensitivity_analysis_params(idx);

    figure(10000 + coeff)
    barh(coeff_values(idx), 'FaceColor', color);

    xlim([-1,1]);
    xticks([-1:0.5:1]);
    xlabel('Sensitivity coefficient');

    ylim([0,length(sensitivity_analysis_params) + 1]);
    yticks([1:1:length(sensitivity_analysis_params)]);
    yticklabels(param_names);
    title(regexp(coeffs(coeff), '_([^_]*)$', 'tokens', 'once'))

    pbaspect([0.75 1 1])

end

%%
toc
%

