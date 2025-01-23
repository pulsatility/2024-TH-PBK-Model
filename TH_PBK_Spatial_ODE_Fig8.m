function dydt = TH_PBK_Spatial_ODE(t, y, option, param)


%% -------------------------- PARAMETERS MAPPING ----------------------------------%%
M = param.M;
L = 100;
dx = L/M;

k1 = param.k1;
k2 = param.k2;
k3 = param.k3;
k4 = param.k4;
k5 = param.k5;
k6 = param.k6;
k7 = param.k7;
k8 = param.k8;
k9 = param.k9;
k10 = param.k10;
k11 = param.k11;
k12 = param.k12;
k20 = param.k20;
k21 = param.k21/M;
k22 = param.k22;
k23 = param.k23/M;
k24 = param.k24;
k25 = param.k25/M;
k26 = param.k26;
k27 = param.k27/M;
k28 = param.k28/M;
k29 = param.k29/M;
k30 = param.k30/M;
k31 = param.k31/M;
k32 = param.k32;
k33 = param.k33;
k34 = param.k34;
k35 = param.k35;
k36 = param.k36;
k37 = param.k37;
k38 = param.k38;
k39 = param.k39;
k42 = param.k42;
k43 = param.k43;
fuT4RBT = param.fuT4RBT;
fuT3RBT = param.fuT3RBT;
fuT4LT  = param.fuT4LT;
fuT3LT  = param.fuT3LT;

TBGtot = param.TBGtot;
TTRtot = param.TTRtot;
ALBtot = param.ALBtot;

QC   = param.QC;
QT   = param.QT;
QRB  = param.QRB;
QL   = param.QL;
VB   = param.VB;
VTB  = param.VTB;
VRBB = param.VRBB;
VRBT = param.VRBT;
VLB  = param.VLB;
VLT  = param.VLT;
%


%% ------------------------- STATE NAME MAPPING----------------------------%%
fT4B    = y(1);
T4TBGB  = y(2);
T4TTRB  = y(3);
T4ALBB  = y(4);
fT3B    = y(5);
T3TBGB  = y(6);
T3TTRB  = y(7);
T3ALBB  = y(8);
TBGB    = y(9);
TTRB    = y(10);
ALBB    = y(11);

fT4T    = y(12);
T4TBGT  = y(13);
T4TTRT  = y(14);
T4ALBT  = y(15);
fT3T    = y(16);
T3TBGT  = y(17);
T3TTRT  = y(18);
T3ALBT  = y(19);
TBGT    = y(20);
TTRT    = y(21);
ALBT    = y(22);

XB      = y(23);
XTTRB   = y(24);
XT      = y(25);
XTTRT   = y(26);
XTBGB   = y(27);
XTBGT   = y(28);

N = 29; %Starting index of the repeating variables
fT4RB   = y(N:N+M-1);
T4TBGRB = y(N+M:N+2*M-1);
T4TTRRB = y(N+2*M:N+3*M-1);
T4ALBRB = y(N+3*M:N+4*M-1);
fT3RB   = y(N+4*M:N+5*M-1);
T3TBGRB = y(N+5*M:N+6*M-1);
T3TTRRB = y(N+6*M:N+7*M-1);
T3ALBRB = y(N+7*M:N+8*M-1);
TBGRB   = y(N+8*M:N+9*M-1);
TTRRB   = y(N+9*M:N+10*M-1);
ALBRB   = y(N+10*M:N+11*M-1);
T4RBT   = y(N+11*M:N+12*M-1);
T3RBT   = y(N+12*M:N+13*M-1);
XRB     = y(N+13*M:N+14*M-1);
XTTRRB  = y(N+14*M:N+15*M-1);
XTBGRB  = y(N+15*M:N+16*M-1);

fT4L    = y(N+16*M:N+17*M-1);
T4TBGL  = y(N+17*M:N+18*M-1);
T4TTRL  = y(N+18*M:N+19*M-1);
T4ALBL  = y(N+19*M:N+20*M-1);
fT3L    = y(N+20*M:N+21*M-1);
T3TBGL  = y(N+21*M:N+22*M-1);
T3TTRL  = y(N+22*M:N+23*M-1);
T3ALBL  = y(N+23*M:N+24*M-1);
TBGL    = y(N+24*M:N+25*M-1);
TTRL    = y(N+25*M:N+26*M-1);
ALBL    = y(N+26*M:N+27*M-1);
T4LT    = y(N+27*M:N+28*M-1);
T3LT    = y(N+28*M:N+29*M-1);
XL      = y(N+29*M:N+30*M-1);
XTTRL   = y(N+30*M:N+31*M-1);
XTBGL   = y(N+31*M:N+32*M-1);
%


%% ------------------------------ ODEs-------------------------------------%%

dydt = zeros(N+32*M-1,1); %make dydt as a column vector as required by MatLab ode function

%Body Blood Compartment 
%fT4B 
dydt(1) = 0; % -k1*fT4B*TBGB + k2*T4TBGB - k3*fT4B*TTRB + k4*T4TTRB - k5*fT4B*ALBB + k6*T4ALBB + (fT4T*QT + fT4RB(M)*QRB + fT4L(M)*QL - fT4B*QC)/VB;

%T4TBGB
dydt(2)	= 0; % k1*fT4B*TBGB - k2*T4TBGB + (T4TBGT*QT + T4TBGRB(M)*QRB + T4TBGL(M)*QL - T4TBGB*QC)/VB;

%T4TTRB
dydt(3) = 0; % k3*fT4B*TTRB - k4*T4TTRB + (T4TTRT*QT + T4TTRRB(M)*QRB + T4TTRL(M)*QL - T4TTRB*QC)/VB;

%T4ALBB
dydt(4) = 0; % k5*fT4B*ALBB - k6*T4ALBB + (T4ALBT*QT + T4ALBRB(M)*QRB + T4ALBL(M)*QL - T4ALBB*QC)/VB;

%fT3B
dydt(5) = 0; % -k7*fT3B*TBGB + k8*T3TBGB - k9*fT3B*TTRB + k10*T3TTRB - k11*fT3B*ALBB + k12*T3ALBB + (fT3T*QT + fT3RB(M)*QRB + fT3L(M)*QL - fT3B*QC)/VB;

%T3TBGB
dydt(6)	= 0; % k7*fT3B*TBGB - k8*T3TBGB + (T3TBGT*QT + T3TBGRB(M)*QRB + T3TBGL(M)*QL - T3TBGB*QC)/VB;

%T3TTRB
dydt(7) = 0; % k9*fT3B*TTRB - k10*T3TTRB + (T3TTRT*QT + T3TTRRB(M)*QRB + T3TTRL(M)*QL - T3TTRB*QC)/VB;

%T3ALBB
dydt(8) = 0; % k11*fT3B*ALBB - k12*T3ALBB + (T3ALBT*QT + T3ALBRB(M)*QRB + T3ALBL(M)*QL - T3ALBB*QC)/VB;

%TBGB
dydt(9) = 0; % -k1*fT4B*TBGB + k2*T4TBGB - k7*fT3B*TBGB + k8*T3TBGB - k42*XB*TBGB + k43*XTBGB + (TBGT*QT + TBGRB(M)*QRB + TBGL(M)*QL - TBGB*QC)/VB;

%TTRB
dydt(10) = 0; % -k3*fT4B*TTRB + k4*T4TTRB - k9*fT3B*TTRB + k10*T3TTRB - k38*XB*TTRB + k39*XTTRB + (TTRT*QT + TTRRB(M)*QRB + TTRL(M)*QL - TTRB*QC)/VB;

%ALBB
dydt(11) = 0; % -k5*fT4B*ALBB + k6*T4ALBB - k11*fT3B*ALBB + k12*T3ALBB + (ALBT*QT + ALBRB(M)*QRB + ALBL(M)*QL - ALBB*QC)/VB;


%Thyroid blood Compartment
%fT4T
dydt(12) = -k1*fT4T*TBGT + k2*T4TBGT - k3*fT4T*TTRT + k4*T4TTRT - k5*fT4T*ALBT + k6*T4ALBT + k20/VTB + (fT4B-fT4T)*QT/VTB;

%T4TBGT
dydt(13) = k1*fT4T*TBGT - k2*T4TBGT + (T4TBGB-T4TBGT)*QT/VTB;

%T4TTRT
dydt(14) = k3*fT4T*TTRT - k4*T4TTRT + (T4TTRB-T4TTRT)*QT/VTB;

%T4ALBT
dydt(15) = k5*fT4T*ALBT - k6*T4ALBT + (T4ALBB-T4ALBT)*QT/VTB;

%fT3T
dydt(16) = -k7*fT3T*TBGT + k8*T3TBGT - k9*fT3T*TTRT + k10*T3TTRT - k11*fT3T*ALBT + k12*T3ALBT + k22/VTB + (fT3B-fT3T)*QT/VTB;

%T3TBGT
dydt(17) = k7*fT3T*TBGT - k8*T3TBGT + (T3TBGB-T3TBGT)*QT/VTB;

%T3TTRT
dydt(18) = k9*fT3T*TTRT - k10*T3TTRT + (T3TTRB-T3TTRT)*QT/VTB;

%T3ALBT
dydt(19) = k11*fT3T*ALBT - k12*T3ALBT + (T3ALBB-T3ALBT)*QT/VTB;

%TBGT
dydt(20) = -k1*fT4T*TBGT + k2*T4TBGT - k7*fT3T*TBGT + k8*T3TBGT - k42*XT*TBGT + k43*XTBGT + (TBGB-TBGT)*QT/VTB;

%TTRT
dydt(21) = -k3*fT4T*TTRT + k4*T4TTRT - k9*fT3T*TTRT + k10*T3TTRT - k38*XT*TTRT + k39*XTTRT + (TTRB-TTRT)*QT/VTB;

%ALBT
dydt(22) = -k5*fT4T*ALBT + k6*T4ALBT - k11*fT3T*ALBT + k12*T3ALBT + (ALBB-ALBT)*QT/VTB;


%Binding Between X and TTR
%XB (Body Blood Compartment)
dydt(23) = k36/VB - k37*XB - k38*XB*TTRB - k42*XB*TBGB + k39*XTTRB + k43*XTBGB + (XT*QT + XRB(M)*QRB + XL(M)*QL - XB*QC)/VB;

%XTTRB (Body Blood Compartment)
dydt(24) = k38*XB*TTRB - k39*XTTRB + (XTTRT*QT + XTTRRB(M)*QRB + XTTRL(M)*QL - XTTRB*QC)/VB;

%XT (Thyroid blood Compartment)
dydt(25) = -k38*XT*TTRT - k42*XT*TBGT + k39*XTTRT + k43*XTBGT + (XB-XT)*QT/VTB;

%XTTRT (Thyroid blood Compartment)
dydt(26) = k38*XT*TTRT - k39*XTTRT + (XTTRB-XTTRT)*QT/VTB;

%XTBGB (Body Blood Compartment)
dydt(27) = k42*XB*TBGB - k43*XTBGB + (XTBGT*QT + XTBGRB(M)*QRB + XTBGL(M)*QL - XTBGB*QC)/VB;

%XTBGT (Thyroid blood Compartment)
dydt(28) = k42*XT*TBGT - k43*XTBGT + (XTBGB-XTBGT)*QT/VTB;

D = 0;  %1E7; %Diffusion constant in tissue blood
E = 0;  %Diffusion constant in tissue propoer

%RB Compartment (first and last segments)
%fT4RB
dydt(N) = D*(1.8E-4)*(fT4RB(2)-fT4RB(1))/(dx^2) - k1*fT4RB(1)*TBGRB(1) + k2*T4TBGRB(1) - k3*fT4RB(1)*TTRRB(1) + k4*T4TTRRB(1) - k5*fT4RB(1)*ALBRB(1) + k6*T4ALBRB(1) + (-k21*fT4RB(1) + k28*T4RBT(1)*fuT4RBT)/(VRBB/M) + (fT4B-fT4RB(1))*QRB/(VRBB/M);
dydt(N+M-1) = D*(1.8E-4)*(fT4RB(M-1)-fT4RB(M))/(dx^2) - k1*fT4RB(M)*TBGRB(M) + k2*T4TBGRB(M) - k3*fT4RB(M)*TTRRB(M) + k4*T4TTRRB(M) - k5*fT4RB(M)*ALBRB(M) + k6*T4ALBRB(M) + (-k21*fT4RB(M) + k28*T4RBT(M)*fuT4RBT)/(VRBB/M) + (fT4RB(M-1)-fT4RB(M))*QRB/(VRBB/M);

%T4TBGRB
dydt(N+M) = D*(1.8E-4)*(T4TBGRB(2)-T4TBGRB(1))/(dx^2) + k1*fT4RB(1)*TBGRB(1) - k2*T4TBGRB(1) + (T4TBGB-T4TBGRB(1))*QRB/(VRBB/M);
dydt(N+2*M-1) = D*(1.8E-4)*(T4TBGRB(M-1)-T4TBGRB(M))/(dx^2) + k1*fT4RB(M)*TBGRB(M) - k2*T4TBGRB(M) + (T4TBGRB(M-1)-T4TBGRB(M))*QRB/(VRBB/M);

%T4TTRRB
dydt(N+2*M) = D*(1.8E-4)*(T4TTRRB(2)-T4TTRRB(1))/(dx^2) + k3*fT4RB(1)*TTRRB(1) - k4*T4TTRRB(1) + (T4TTRB-T4TTRRB(1))*QRB/(VRBB/M);
dydt(N+3*M-1) = D*(1.8E-4)*(T4TTRRB(M-1)-T4TTRRB(M))/(dx^2) + k3*fT4RB(M)*TTRRB(M) - k4*T4TTRRB(M) + (T4TTRRB(M-1)-T4TTRRB(M))*QRB/(VRBB/M);

%T4ALBRB
dydt(N+3*M) = D*(1.8E-4)*(T4ALBRB(2)-T4ALBRB(1))/(dx^2) + k5*fT4RB(1)*ALBRB(1) - k6*T4ALBRB(1) + (T4ALBB-T4ALBRB(1))*QRB/(VRBB/M);
dydt(N+4*M-1) = D*(1.8E-4)*(T4ALBRB(M-1)-T4ALBRB(M))/(dx^2) + k5*fT4RB(M)*ALBRB(M) - k6*T4ALBRB(M) + (T4ALBRB(M-1)-T4ALBRB(M))*QRB/(VRBB/M);

%fT3RB
dydt(N+4*M) = D*(1.8E-4)*(fT3RB(2)-fT3RB(1))/(dx^2) - k7*fT3RB(1)*TBGRB(1) + k8*T3TBGRB(1) - k9*fT3RB(1)*TTRRB(1) + k10*T3TTRRB(1) - k11*fT3RB(1)*ALBRB(1) + k12*T3ALBRB(1) + (-k23*fT3RB(1) + k29*T3RBT(1)*fuT3RBT)/(VRBB/M) + (fT3B-fT3RB(1))*QRB/(VRBB/M);
dydt(N+5*M-1) = D*(1.8E-4)*(fT3RB(M-1)-fT3RB(M))/(dx^2) - k7*fT3RB(M)*TBGRB(M) + k8*T3TBGRB(M) - k9*fT3RB(M)*TTRRB(M) + k10*T3TTRRB(M) - k11*fT3RB(M)*ALBRB(M) + k12*T3ALBRB(M) + (-k23*fT3RB(M) + k29*T3RBT(M)*fuT3RBT)/(VRBB/M) + (fT3RB(M-1)-fT3RB(M))*QRB/(VRBB/M);

%T3TBGRB
dydt(N+5*M) = D*(1.8E-4)*(T3TBGRB(2)-T3TBGRB(1))/(dx^2) + k7*fT3RB(1)*TBGRB(1) - k8*T3TBGRB(1) + (T3TBGB-T3TBGRB(1))*QRB/(VRBB/M);
dydt(N+6*M-1) = D*(1.8E-4)*(T3TBGRB(M-1)-T3TBGRB(M))/(dx^2) + k7*fT3RB(M)*TBGRB(M) - k8*T3TBGRB(M) + (T3TBGRB(M-1)-T3TBGRB(M))*QRB/(VRBB/M);

%T3TTRRB
dydt(N+6*M) = D*(1.8E-4)*(T3TTRRB(2)-T3TTRRB(1))/(dx^2) + k9*fT3RB(1)*TTRRB(1) - k10*T3TTRRB(1) + (T3TTRB-T3TTRRB(1))*QRB/(VRBB/M);
dydt(N+7*M-1) = D*(1.8E-4)*(T3TTRRB(M-1)-T3TTRRB(M))/(dx^2) + k9*fT3RB(M)*TTRRB(M) - k10*T3TTRRB(M) + (T3TTRRB(M-1)-T3TTRRB(M))*QRB/(VRBB/M);

%T3ALBRB
dydt(N+7*M) = D*(1.8E-4)*(T3ALBRB(2)-T3ALBRB(1))/(dx^2) + k11*fT3RB(1)*ALBRB(1) - k12*T3ALBRB(1) + (T3ALBB-T3ALBRB(1))*QRB/(VRBB/M);
dydt(N+8*M-1) = D*(1.8E-4)*(T3ALBRB(M-1)-T3ALBRB(M))/(dx^2) + k11*fT3RB(M)*ALBRB(M) - k12*T3ALBRB(M) + (T3ALBRB(M-1)-T3ALBRB(M))*QRB/(VRBB/M);

%TBGRB
dydt(N+8*M) = D*(1.8E-4)*(TBGRB(2)-TBGRB(1))/(dx^2) - k1*fT4RB(1)*TBGRB(1) + k2*T4TBGRB(1) - k7*fT3RB(1)*TBGRB(1) + k8*T3TBGRB(1) - k42*XRB(1)*TBGRB(1) + k43*XTBGRB(1) + (TBGB-TBGRB(1))*QRB/(VRBB/M);
dydt(N+9*M-1) = D*(1.8E-4)*(TBGRB(M-1)-TBGRB(M))/(dx^2) - k1*fT4RB(M)*TBGRB(M) + k2*T4TBGRB(M) - k7*fT3RB(M)*TBGRB(M) + k8*T3TBGRB(M) - k42*XRB(M)*TBGRB(M) + k43*XTBGRB(M) + (TBGRB(M-1)-TBGRB(M))*QRB/(VRBB/M);

%TTRRB
dydt(N+9*M) = D*(1.8E-4)*(TTRRB(2)-TTRRB(1))/(dx^2) - k3*fT4RB(1)*TTRRB(1) + k4*T4TTRRB(1) - k9*fT3RB(1)*TTRRB(1) + k10*T3TTRRB(1) - k38*XRB(1)*TTRRB(1) + k39*XTTRRB(1) + (TTRB-TTRRB(1))*QRB/(VRBB/M);
dydt(N+10*M-1) = D*(1.8E-4)*(TTRRB(M-1)-TTRRB(M))/(dx^2) - k3*fT4RB(M)*TTRRB(M) + k4*T4TTRRB(M) - k9*fT3RB(M)*TTRRB(M) + k10*T3TTRRB(M) - k38*XRB(M)*TTRRB(M) + k39*XTTRRB(M) + (TTRRB(M-1)-TTRRB(M))*QRB/(VRBB/M);

%ALBRB
dydt(N+10*M) = D*(1.8E-4)*(ALBRB(2)-ALBRB(1))/(dx^2) - k5*fT4RB(1)*ALBRB(1) + k6*T4ALBRB(1) - k11*fT3RB(1)*ALBRB(1) + k12*T3ALBRB(1) + (ALBB-ALBRB(1))*QRB/(VRBB/M);
dydt(N+11*M-1) = D*(1.8E-4)*(ALBRB(M-1)-ALBRB(M))/(dx^2) - k5*fT4RB(M)*ALBRB(M) + k6*T4ALBRB(M) - k11*fT3RB(M)*ALBRB(M) + k12*T3ALBRB(M) + (ALBRB(M-1)-ALBRB(M))*QRB/(VRBB/M);

%T4RBT
dydt(N+11*M) = E*(1.8E-4)*(T4RBT(2)-T4RBT(1))/(dx^2) + (k21*fT4RB(1) - k28*T4RBT(1)*fuT4RBT)/(VRBT/M) - k24*T4RBT(1)*fuT4RBT - k32*T4RBT(1)*fuT4RBT;
dydt(N+12*M-1) = E*(1.8E-4)*(T4RBT(M-1)-T4RBT(M))/(dx^2) + (k21*fT4RB(M) - k28*T4RBT(M)*fuT4RBT)/(VRBT/M) - k24*T4RBT(M)*fuT4RBT - k32*T4RBT(M)*fuT4RBT;

%T3RBT
dydt(N+12*M) = E*(1.8E-4)*(T3RBT(2)-T3RBT(1))/(dx^2) + (k23*fT3RB(1) - k29*T3RBT(1)*fuT3RBT)/(VRBT/M) + k24*T4RBT(1)*fuT4RBT - k33*T3RBT(1)*fuT3RBT;
dydt(N+13*M-1) = E*(1.8E-4)*(T3RBT(M-1)-T3RBT(M))/(dx^2) + (k23*fT3RB(M) - k29*T3RBT(M)*fuT3RBT)/(VRBT/M) + k24*T4RBT(M)*fuT4RBT - k33*T3RBT(M)*fuT3RBT;

%XRB
dydt(N+13*M) = D*(1.8E-4)*(XRB(2)-XRB(1))/(dx^2) - k38*XRB(1)*TTRRB(1) - k42*XRB(1)*TBGRB(1) + k39*XTTRRB(1) + k43*XTBGRB(1) + (XB-XRB(1))*QRB/(VRBB/M);
dydt(N+14*M-1) = D*(1.8E-4)*(XRB(M-1)-XRB(M))/(dx^2) - k38*XRB(M)*TTRRB(M) - k42*XRB(M)*TBGRB(M) + k39*XTTRRB(M) + k43*XTBGRB(M) + (XRB(M-1)-XRB(M))*QRB/(VRBB/M);

%XTTRRB
dydt(N+14*M) = D*(1.8E-4)*(XTTRRB(2)-XTTRRB(1))/(dx^2) + k38*XRB(1)*TTRRB(1) - k39*XTTRRB(1) + (XTTRB-XTTRRB(1))*QRB/(VRBB/M);
dydt(N+15*M-1) = D*(1.8E-4)*(XTTRRB(M-1)-XTTRRB(M))/(dx^2) + k38*XRB(M)*TTRRB(M) - k39*XTTRRB(M) + (XTTRRB(M-1)-XTTRRB(M))*QRB/(VRBB/M);

%XTBGRB
dydt(N+15*M) = D*(1.8E-4)*(XTBGRB(2)-XTBGRB(1))/(dx^2) + k42*XRB(1)*TBGRB(1) - k43*XTBGRB(1) + (XTBGB-XTBGRB(1))*QRB/(VRBB/M);
dydt(N+16*M-1) = D*(1.8E-4)*(XTBGRB(M-1)-XTBGRB(M))/(dx^2) + k42*XRB(M)*TBGRB(M) - k43*XTBGRB(M) + (XTBGRB(M-1)-XTBGRB(M))*QRB/(VRBB/M);


%Liver Compartment (first and last segments)
%fT4L
dydt(N+16*M) = D*(1.8E-4)*(fT4L(2)-fT4L(1))/(dx^2) - k1*fT4L(1)*TBGL(1) + k2*T4TBGL(1) - k3*fT4L(1)*TTRL(1) + k4*T4TTRL(1) - k5*fT4L(1)*ALBL(1) + k6*T4ALBL(1) + (-k25*fT4L(1) + k30*T4LT(1)*fuT4LT)/(VLB/M) + (fT4B-fT4L(1))*QL/(VLB/M);
dydt(N+17*M-1) = D*(1.8E-4)*(fT4L(M-1)-fT4L(M))/(dx^2) - k1*fT4L(M)*TBGL(M) + k2*T4TBGL(M) - k3*fT4L(M)*TTRL(M) + k4*T4TTRL(M) - k5*fT4L(M)*ALBL(M) + k6*T4ALBL(M) + (-k25*fT4L(M) + k30*T4LT(M)*fuT4LT)/(VLB/M) + (fT4L(M-1)-fT4L(M))*QL/(VLB/M);

%T4TBGL
dydt(N+17*M) = D*(1.8E-4)*(T4TBGL(2)-T4TBGL(1))/(dx^2) + k1*fT4L(1)*TBGL(1) - k2*T4TBGL(1) + (T4TBGB-T4TBGL(1))*QL/(VLB/M);
dydt(N+18*M-1) = D*(1.8E-4)*(T4TBGL(M-1)-T4TBGL(M))/(dx^2) + k1*fT4L(M)*TBGL(M) - k2*T4TBGL(M) + (T4TBGL(M-1)-T4TBGL(M))*QL/(VLB/M);

%T4TTRL
dydt(N+18*M) = D*(1.8E-4)*(T4TTRL(2)-T4TTRL(1))/(dx^2) + k3*fT4L(1)*TTRL(1) - k4*T4TTRL(1) + (T4TTRB-T4TTRL(1))*QL/(VLB/M);
dydt(N+19*M-1) = D*(1.8E-4)*(T4TTRL(M-1)-T4TTRL(M))/(dx^2) + k3*fT4L(M)*TTRL(M) - k4*T4TTRL(M) + (T4TTRL(M-1)-T4TTRL(M))*QL/(VLB/M);

%T4ALBL
dydt(N+19*M) = D*(1.8E-4)*(T4ALBL(2)-T4ALBL(1))/(dx^2) + k5*fT4L(1)*ALBL(1) - k6*T4ALBL(1) + (T4ALBB-T4ALBL(1))*QL/(VLB/M);
dydt(N+20*M-1) = D*(1.8E-4)*(T4ALBL(M-1)-T4ALBL(M))/(dx^2) + k5*fT4L(M)*ALBL(M) - k6*T4ALBL(M) + (T4ALBL(M-1)-T4ALBL(M))*QL/(VLB/M);

%fT3L
dydt(N+20*M) = D*(1.8E-4)*(fT3L(2)-fT3L(1))/(dx^2) - k7*fT3L(1)*TBGL(1) + k8*T3TBGL(1) - k9*fT3L(1)*TTRL(1) + k10*T3TTRL(1) - k11*fT3L(1)*ALBL(1) + k12*T3ALBL(1) + (-k27*fT3L(1) + k31*T3LT(1)*fuT3LT)/(VLB/M) + (fT3B-fT3L(1))*QL/(VLB/M);
dydt(N+21*M-1) = D*(1.8E-4)*(fT3L(M-1)-fT3L(M))/(dx^2) - k7*fT3L(M)*TBGL(M) + k8*T3TBGL(M) - k9*fT3L(M)*TTRL(M) + k10*T3TTRL(M) - k11*fT3L(M)*ALBL(M) + k12*T3ALBL(M) + (-k27*fT3L(M) + k31*T3LT(M)*fuT3LT)/(VLB/M) + (fT3L(M-1)-fT3L(M))*QL/(VLB/M);

%T3TBGL
dydt(N+21*M) = D*(1.8E-4)*(T3TBGL(2)-T3TBGL(1))/(dx^2) + k7*fT3L(1)*TBGL(1) - k8*T3TBGL(1) + (T3TBGB-T3TBGL(1))*QL/(VLB/M);
dydt(N+22*M-1) = D*(1.8E-4)*(T3TBGL(M-1)-T3TBGL(M))/(dx^2) + k7*fT3L(M)*TBGL(M) - k8*T3TBGL(M) + (T3TBGL(M-1)-T3TBGL(M))*QL/(VLB/M);

%T3TTRL
dydt(N+22*M) = D*(1.8E-4)*(T3TTRL(2)-T3TTRL(1))/(dx^2) + k9*fT3L(1)*TTRL(1) - k10*T3TTRL(1) + (T3TTRB-T3TTRL(1))*QL/(VLB/M);
dydt(N+23*M-1) = D*(1.8E-4)*(T3TTRL(M-1)-T3TTRL(M))/(dx^2) + k9*fT3L(M)*TTRL(M) - k10*T3TTRL(M) + (T3TTRL(M-1)-T3TTRL(M))*QL/(VLB/M);

%T3ALBL
dydt(N+23*M) = D*(1.8E-4)*(T3ALBL(2)-T3ALBL(1))/(dx^2) + k11*fT3L(1)*ALBL(1) - k12*T3ALBL(1) + (T3ALBB-T3ALBL(1))*QL/(VLB/M);
dydt(N+24*M-1) = D*(1.8E-4)*(T3ALBL(M-1)-T3ALBL(M))/(dx^2) + k11*fT3L(M)*ALBL(M) - k12*T3ALBL(M) + (T3ALBL(M-1)-T3ALBL(M))*QL/(VLB/M);

%TBGL
dydt(N+24*M) = D*(1.8E-4)*(TBGL(2)-TBGL(1))/(dx^2) - k1*fT4L(1)*TBGL(1) + k2*T4TBGL(1) - k7*fT3L(1)*TBGL(1) + k8*T3TBGL(1) - k42*XL(1)*TBGL(1) + k43*XTBGL(1) + (TBGB-TBGL(1))*QL/(VLB/M);
dydt(N+25*M-1) = D*(1.8E-4)*(TBGL(M-1)-TBGL(M))/(dx^2) - k1*fT4L(M)*TBGL(M) + k2*T4TBGL(M) - k7*fT3L(M)*TBGL(M) + k8*T3TBGL(M) - k42*XL(M)*TBGL(M) + k43*XTBGL(M) + (TBGL(M-1)-TBGL(M))*QL/(VLB/M);

%TTRL
dydt(N+25*M) = D*(1.8E-4)*(TTRRB(2)-TTRRB(1))/(dx^2) - k3*fT4L(1)*TTRL(1) + k4*T4TTRL(1) - k9*fT3L(1)*TTRL(1) + k10*T3TTRL(1) - k38*XL(1)*TTRL(1) + k39*XTTRL(1) + (TTRB-TTRL(1))*QL/(VLB/M);
dydt(N+26*M-1) = D*(1.8E-4)*(TTRRB(M-1)-TTRRB(M))/(dx^2) - k3*fT4L(M)*TTRL(M) + k4*T4TTRL(M) - k9*fT3L(M)*TTRL(M) + k10*T3TTRL(M) - k38*XL(M)*TTRL(M) + k39*XTTRL(M) + (TTRL(M-1)-TTRL(M))*QL/(VLB/M);

%ALBL
dydt(N+26*M) = D*(1.8E-4)*(ALBL(2)-ALBL(1))/(dx^2) - k5*fT4L(1)*ALBL(1) + k6*T4ALBL(1) - k11*fT3L(1)*ALBL(1) + k12*T3ALBL(1) + (ALBB-ALBL(1))*QL/(VLB/M);
dydt(N+27*M-1) = D*(1.8E-4)*(ALBL(M-1)-ALBL(M))/(dx^2) - k5*fT4L(M)*ALBL(M) + k6*T4ALBL(M) - k11*fT3L(M)*ALBL(M) + k12*T3ALBL(M) + (ALBL(M-1)-ALBL(M))*QL/(VLB/M);

%T4LT
dydt(N+27*M) = E*(1.8E-4)*(T4LT(2)-T4LT(1))/(dx^2) + (k25*fT4L(1) - k30*T4LT(1)*fuT4LT)/(VLT/M) - k26*T4LT(1)*fuT4LT - k34*T4LT(1)*fuT4LT;
dydt(N+28*M-1) = E*(1.8E-4)*(T4LT(M-1)-T4LT(M))/(dx^2) + (k25*fT4L(M) - k30*T4LT(M)*fuT4LT)/(VLT/M) - k26*T4LT(M)*fuT4LT - k34*T4LT(M)*fuT4LT;

%T3LT
dydt(N+28*M) = E*(1.8E-4)*(T3LT(2)-T3LT(1))/(dx^2) + (k27*fT3L(1) - k31*T3LT(1)*fuT3LT)/(VLT/M) + k26*T4LT(1)*fuT4LT - k35*T3LT(1)*fuT3LT;
dydt(N+29*M-1) = E*(1.8E-4)*(T3LT(M-1)-T3LT(M))/(dx^2) + (k27*fT3L(M) - k31*T3LT(M)*fuT3LT)/(VLT/M) + k26*T4LT(M)*fuT4LT - k35*T3LT(M)*fuT3LT;

%XL
dydt(N+29*M) = D*(1.8E-4)*(XL(2)-XL(1))/(dx^2) - k38*XL(1)*TTRL(1) - k42*XL(1)*TBGL(1) + k39*XTTRL(1) + k43*XTBGL(1) + (XB-XL(1))*QL/(VLB/M);
dydt(N+30*M-1) = D*(1.8E-4)*(XL(M-1)-XL(M))/(dx^2) - k38*XL(M)*TTRL(M) - k42*XL(M)*TBGL(M) + k39*XTTRL(M) + k43*XTBGL(M) + (XL(M-1)-XL(M))*QL/(VLB/M);

%XTTRL
dydt(N+30*M) = D*(1.8E-4)*(XTTRL(2)-XTTRL(1))/(dx^2) + k38*XL(1)*TTRL(1) - k39*XTTRL(1) + (XTTRB-XTTRL(1))*QL/(VLB/M);
dydt(N+31*M-1) = D*(1.8E-4)*(XTTRL(M-1)-XTTRL(M))/(dx^2) + k38*XL(M)*TTRL(M) - k39*XTTRL(M) + (XTTRL(M-1)-XTTRL(M))*QL/(VLB/M);

%XTBGL
dydt(N+31*M) = D*(1.8E-4)*(XTBGL(2)-XTBGL(1))/(dx^2) + k42*XL(1)*TBGL(1) - k43*XTBGL(1) + (XTBGB-XTBGL(1))*QL/(VLB/M);
dydt(N+32*M-1) = D*(1.8E-4)*(XTBGL(M-1)-XTBGL(M))/(dx^2) + k42*XL(M)*TBGL(M) - k43*XTBGL(M) + (XTBGL(M-1)-XTBGL(M))*QL/(VLB/M);

%Intermediate segments for RB and Liver Compartments
for i=2:(M-1)
    dydt(i-1+N) = D*(1.8E-4)*(fT4RB(i+1)+fT4RB(i-1)-2*fT4RB(i))/(dx^2) - k1*fT4RB(i)*TBGRB(i) + k2*T4TBGRB(i) - k3*fT4RB(i)*TTRRB(i) + k4*T4TTRRB(i) - k5*fT4RB(i)*ALBRB(i) + k6*T4ALBRB(i) + (-k21*fT4RB(i) + k28*T4RBT(i)*fuT4RBT)/(VRBB/M) + ((fT4RB(i-1)-fT4RB(i))*QRB)/(VRBB/M);
    dydt(i+M-1+N) = D*(1.8E-4)*(T4TBGRB(i+1)+T4TBGRB(i-1)-2*T4TBGRB(i))/(dx^2) + k1*fT4RB(i)*TBGRB(i) - k2*T4TBGRB(i) + (T4TBGRB(i-1)-T4TBGRB(i))*QRB/(VRBB/M);
    dydt(i+2*M-1+N) = D*(1.8E-4)*(T4TTRRB(i+1)+T4TTRRB(i-1)-2*T4TTRRB(i))/(dx^2) + k3*fT4RB(i)*TTRRB(i) - k4*T4TTRRB(i) + (T4TTRRB(i-1)-T4TTRRB(i))*QRB/(VRBB/M);
    dydt(i+3*M-1+N) = D*(1.8E-4)*(T4ALBRB(i+1)+T4ALBRB(i-1)-2*T4ALBRB(i))/(dx^2) + k5*fT4RB(i)*ALBRB(i) - k6*T4ALBRB(i) + (T4ALBRB(i-1)-T4ALBRB(i))*QRB/(VRBB/M);
    dydt(i+4*M-1+N) = D*(1.8E-4)*(fT3RB(i+1)+fT3RB(i-1)-2*fT3RB(i))/(dx^2) - k7*fT3RB(i)*TBGRB(i) + k8*T3TBGRB(i) - k9*fT3RB(i)*TTRRB(i) + k10*T3TTRRB(i) - k11*fT3RB(i)*ALBRB(i) + k12*T3ALBRB(i) + (-k23*fT3RB(i) + k29*T3RBT(i)*fuT3RBT)/(VRBB/M) + (fT3RB(i-1)-fT3RB(i))*QRB/(VRBB/M);
    dydt(i+5*M-1+N) = D*(1.8E-4)*(T3TBGRB(i+1)+T3TBGRB(i-1)-2*T3TBGRB(i))/(dx^2) + k7*fT3RB(i)*TBGRB(i) - k8*T3TBGRB(i) + (T3TBGRB(i-1)-T3TBGRB(i))*QRB/(VRBB/M);
    dydt(i+6*M-1+N) = D*(1.8E-4)*(T3TTRRB(i+1)+T3TTRRB(i-1)-2*T3TTRRB(i))/(dx^2) + k9*fT3RB(i)*TTRRB(i) - k10*T3TTRRB(i) + (T3TTRRB(i-1)-T3TTRRB(i))*QRB/(VRBB/M);
    dydt(i+7*M-1+N) = D*(1.8E-4)*(T3ALBRB(i+1)+T3ALBRB(i-1)-2*T3ALBRB(i))/(dx^2) + k11*fT3RB(i)*ALBRB(i) - k12*T3ALBRB(i) + (T3ALBRB(i-1)-T3ALBRB(i))*QRB/(VRBB/M); 
    dydt(i+8*M-1+N) = D*(1.8E-4)*(TBGRB(i+1)+TBGRB(i-1)-2*TBGRB(i))/(dx^2) - k1*fT4RB(i)*TBGRB(i) + k2*T4TBGRB(i) - k7*fT3RB(i)*TBGRB(i) + k8*T3TBGRB(i) - k42*XRB(i)*TBGRB(i) + k43*XTBGRB(i) + (TBGRB(i-1)-TBGRB(i))*QRB/(VRBB/M);
    dydt(i+9*M-1+N) = D*(1.8E-4)*(TTRRB(i+1)+TTRRB(i-1)-2*TTRRB(i))/(dx^2) - k3*fT4RB(i)*TTRRB(i) + k4*T4TTRRB(i) - k9*fT3RB(i)*TTRRB(i) + k10*T3TTRRB(i) - k38*XRB(i)*TTRRB(i) + k39*XTTRRB(i) + (TTRRB(i-1)-TTRRB(i))*QRB/(VRBB/M);
    dydt(i+10*M-1+N) = D*(1.8E-4)*(ALBRB(i+1)+ALBRB(i-1)-2*ALBRB(i))/(dx^2) - k5*fT4RB(i)*ALBRB(i) + k6*T4ALBRB(i) - k11*fT3RB(i)*ALBRB(i) + k12*T3ALBRB(i) + (ALBRB(i-1)-ALBRB(i))*QRB/(VRBB/M);
    dydt(i+11*M-1+N) = E*(1.8E-4)*(T4RBT(i+1)+T4RBT(i-1)-2*T4RBT(i))/(dx^2) + (k21*fT4RB(i) - k28*T4RBT(i)*fuT4RBT)/(VRBT/M) - k24*T4RBT(i)*fuT4RBT - k32*T4RBT(i)*fuT4RBT;
    dydt(i+12*M-1+N) = E*(1.8E-4)*(T3RBT(i+1)+T3RBT(i-1)-2*T3RBT(i))/(dx^2) + (k23*fT3RB(i) - k29*T3RBT(i)*fuT3RBT)/(VRBT/M) + k24*T4RBT(i)*fuT4RBT - k33*T3RBT(i)*fuT3RBT;
    dydt(i+13*M-1+N) = D*(1.8E-4)*(XRB(i+1)+XRB(i-1)-2*XRB(i))/(dx^2) - k38*XRB(i)*TTRRB(i) + k39*XTTRRB(i) + (XRB(i-1)-XRB(i))*QRB/(VRBB/M);
    dydt(i+14*M-1+N) = D*(1.8E-4)*(XTTRRB(i+1)+XTTRRB(i-1)-2*XTTRRB(i))/(dx^2) + k38*XRB(i)*TTRRB(i) - k39*XTTRRB(i) + (XTTRRB(i-1)-XTTRRB(i))*QRB/(VRBB/M);
    dydt(i+15*M-1+N) = D*(1.8E-4)*(XTBGRB(i+1)+XTBGRB(i-1)-2*XTBGRB(i))/(dx^2) + k42*XRB(i)*TBGRB(i) - k43*XTBGRB(i) + (XTBGRB(i-1)-XTBGRB(i))*QRB/(VRBB/M);
    dydt(i+16*M-1+N) = D*(1.8E-4)*(fT4L(i+1)+fT4L(i-1)-2*fT4L(i))/(dx^2) - k1*fT4L(i)*TBGL(i) + k2*T4TBGL(i) - k3*fT4L(i)*TTRL(i) + k4*T4TTRL(i) - k5*fT4L(i)*ALBL(i) + k6*T4ALBL(i) + (-k25*fT4L(i) + k30*T4LT(i)*fuT4LT)/(VLB/M) + (fT4L(i-1)-fT4L(i))*QL/(VLB/M);
    dydt(i+17*M-1+N) = D*(1.8E-4)*(T4TBGL(i+1)+T4TBGL(i-1)-2*T4TBGL(i))/(dx^2) + k1*fT4L(i)*TBGL(i) - k2*T4TBGL(i) + (T4TBGL(i-1)-T4TBGL(i))*QL/(VLB/M);
    dydt(i+18*M-1+N) = D*(1.8E-4)*(T4TTRL(i+1)+T4TTRL(i-1)-2*T4TTRL(i))/(dx^2) + k3*fT4L(i)*TTRL(i) - k4*T4TTRL(i) + (T4TTRL(i-1)-T4TTRL(i))*QL/(VLB/M);
    dydt(i+19*M-1+N) = D*(1.8E-4)*(T4ALBL(i+1)+T4ALBL(i-1)-2*T4ALBL(i))/(dx^2) + k5*fT4L(i)*ALBL(i) - k6*T4ALBL(i) + (T4ALBL(i-1)-T4ALBL(i))*QL/(VLB/M);
    dydt(i+20*M-1+N) = D*(1.8E-4)*(fT3L(i+1)+fT3L(i-1)-2*fT3L(i))/(dx^2) - k7*fT3L(i)*TBGL(i) + k8*T3TBGL(i) - k9*fT3L(i)*TTRL(i) + k10*T3TTRL(i) - k11*fT3L(i)*ALBL(i) + k12*T3ALBL(i) + (-k27*fT3L(i) + k31*T3LT(i)*fuT3LT)/(VLB/M) + (fT3L(i-1)-fT3L(i))*QL/(VLB/M);
    dydt(i+21*M-1+N) = D*(1.8E-4)*(T3TBGL(i+1)+T3TBGL(i-1)-2*T3TBGL(i))/(dx^2) + k7*fT3L(i)*TBGL(i) - k8*T3TBGL(i) + (T3TBGL(i-1)-T3TBGL(i))*QL/(VLB/M);
    dydt(i+22*M-1+N) = D*(1.8E-4)*(T3TTRL(i+1)+T3TTRL(i-1)-2*T3TTRL(i))/(dx^2) + k9*fT3L(i)*TTRL(i) - k10*T3TTRL(i) + (T3TTRL(i-1)-T3TTRL(i))*QL/(VLB/M);
    dydt(i+23*M-1+N) = D*(1.8E-4)*(T3ALBL(i+1)+T3ALBL(i-1)-2*T3ALBL(i))/(dx^2) + k11*fT3L(i)*ALBL(i) - k12*T3ALBL(i) + (T3ALBL(i-1)-T3ALBL(i))*QL/(VLB/M);
    dydt(i+24*M-1+N) = D*(1.8E-4)*(TBGL(i+1)+TBGL(i-1)-2*TBGL(i))/(dx^2) - k1*fT4L(i)*TBGL(i) + k2*T4TBGL(i) - k7*fT3L(i)*TBGL(i) + k8*T3TBGL(i) - k42*XL(i)*TBGL(i) + k43*XTBGL(i) + (TBGL(i-1)-TBGL(i))*QL/(VLB/M);
    dydt(i+25*M-1+N) = D*(1.8E-4)*(TTRL(i+1)+TTRL(i-1)-2*TTRL(i))/(dx^2) - k3*fT4L(i)*TTRL(i) + k4*T4TTRL(i) - k9*fT3L(i)*TTRL(i) + k10*T3TTRL(i) - k38*XL(i)*TTRL(i) + k39*XTTRL(i) + (TTRL(i-1)-TTRL(i))*QL/(VLB/M);
    dydt(i+26*M-1+N) = D*(1.8E-4)*(ALBL(i+1)+ALBL(i-1)-2*ALBL(i))/(dx^2) - k5*fT4L(i)*ALBL(i) + k6*T4ALBL(i) - k11*fT3L(i)*ALBL(i) + k12*T3ALBL(i) + (ALBL(i-1)-ALBL(i))*QL/(VLB/M);
    dydt(i+27*M-1+N) = E*(1.8E-4)*(T4LT(i+1)+T4LT(i-1)-2*T4LT(i))/(dx^2) + (k25*fT4L(i) - k30*T4LT(i)*fuT4LT)/(VLT/M) - k26*T4LT(i)*fuT4LT - k34*T4LT(i)*fuT4LT;
    dydt(i+28*M-1+N) = E*(1.8E-4)*(T3LT(i+1)+T3LT(i-1)-2*T3LT(i))/(dx^2) + (k27*fT3L(i) - k31*T3LT(i)*fuT3LT)/(VLT/M) + k26*T4LT(i)*fuT4LT - k35*T3LT(i)*fuT3LT;
    dydt(i+29*M-1+N) = D*(1.8E-4)*(XL(i+1)+XL(i-1)-2*XL(i))/(dx^2) - k38*XL(i)*TTRL(i) - k42*XL(i)*TBGL(i) + k39*XTTRL(i) + k43*XTBGL(i) + (XL(i-1)-XL(i))*QL/(VLB/M);
    dydt(i+30*M-1+N) = D*(1.8E-4)*(XTTRL(i+1)+XTTRL(i-1)-2*XTTRL(i))/(dx^2) + k38*XL(i)*TTRL(i) - k39*XTTRL(i) + (XTTRL(i-1)-XTTRL(i))*QL/(VLB/M);
    dydt(i+31*M-1+N) = D*(1.8E-4)*(XTBGL(i+1)+XTBGL(i-1)-2*XTBGL(i))/(dx^2) + k42*XL(i)*TBGL(i) - k43*XTBGL(i) + (XTBGL(i-1)-XTBGL(i))*QL/(VLB/M);
end
  
end