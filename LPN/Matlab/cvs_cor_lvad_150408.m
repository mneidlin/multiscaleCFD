function [RES CO]=cvs_cor_lvad_150408(a,b,M)
                                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                              %     Michael Neidlin       %
                                              %        08-04-2015         %
                                              % michael.neidlin@gmail.com %
                                              %                           %
                                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Euler solution of entire LP CVS (coronary circulation with and without LVAD), as seen in Simulink file!
%Input 
% - a and b are the left and right endpoints
% - ya is the initial condition y(a)
% - M is the number of steps
%Output - E=[T' Y'] where T is the vector of abscissas and
% Y is the vector of ordinates
% Forward Euler Explicit: y(t+1)=y(t)+h*ODE with ODE=> dy/dt=-y-3*t

% CO is mean systemic output (Heart + VAD)
% RES is result matrix defined in the last row of function

%%%%%%%%%%%%%%%%%%%%+
%TIMESTEP OF 1e-4 is necessary for Explicit! TS of 2e-3 up to 5e-3 is possible for implicit/explicit mix !!! 
%values of R (mmHg ml/s) and C ml/mmHg and CQi (ml /s mmHg^0.5) and Volumes
%of Compartments V0 (ml)
%1: sart_UB 2: svn_UB 3: sart_lb, 4:svn_LB, 5:RA, 6:RV, 7:pulm_artery 8:
%part 9: pvn 10: LA 11:LV 12:AorticArch 13:LCORat 14:LCOR 15:LVAD
%%%%%%%%%%%%%%%%%%%%
%PARAMETERS ONLY FOR ONE COMPARTMENT OF CORONARIES
%NONLINEAR AORTIC VALVE IS INCLUDED
%Parameters were manually optimised to fit physiological values (/x *x)
RsartUB=2.96;
RsvnUB=0.33;
RsartLB=1.44;
RsvnLB=0.05;
CQtri=400;
CQpa=350;
Rpa=0.02; 
Rpart=0.31/2;
Rpvn=0.075;
CQmi=400;
CQao=350;
Rao=0.07; %THIS IS THE SUM OF PARALLEL RESISTANCES Zc1,Zc2 from Multiscale model!!!
Rlcor=4.5;
Rcan_in=0.05;
Rcan_out=0.05;


%Valve parameters
A_real=0;
Kav=5000;
Kp=50;
Angle_max=85;

CsartUB=1.35*0.85; %1.35*0.8 for longterm stable simulation
CsvnUB=4;
CsartLB=1.95*0.6; %1.95*0.5 for longterm stable simulation
CsvnLB=15;
EmaxRA=0.225;
EdRA=0.1575;
EmaxRV=0.45;
Cpa=0.18;
Cpart=3.8;
Cpvn=20.5; 
EmaxLA=0.45;
EdLA=0.2125;
EmaxLV=2.8; %CHANGE 2.8=>0.6 LHF
Cao=0.3;
Clcor=0.05;
Crcor=0;




V0sartUB=80;
V0svnUB=60;
V0sartLB=100;
V0svnLB=150;%CHANGE 150 => 80 for LHF
V0RA=20;%40;
V0RV=145; 
V0pa=10;%10;
V0part=155;%155;
V0pvn=260;%260; %CHANGE 260 => 320 for LHF
V0LA=22;
V0LV=105; 
V0ao=15;
V0lcor=0.1;
V0rcor=0.1;


R(1)=RsartUB;
R(2)=RsvnUB;
R(3)=RsartLB;
R(4)=RsvnLB;
R(7)=Rpa;
R(8)=Rpart;
R(9)=Rpvn;
R(12)=Rao;
R(13)=Rlcor;

R(15)=Rcan_in;
R(16)=Rcan_out;

C(1)=CsartUB;
C(2)=CsvnUB;
C(3)=CsartLB;
C(4)=CsvnLB;
C(7)=Cpa;
C(8)=Cpart;
C(9)=Cpvn;
C(12)=Cao;
C(13)=Clcor;
C(14)=Crcor;

LV_Pd_beta=3;%2.5;        %Constants for LVEDPVR
LV_Pd_kappa=0.0396;%0.033;
LV_Pd_alpha=0.0768;%0.064;

RV_Pd_beta=3;%2.5;        %Constants for RVEDPVR
RV_Pd_kappa=0.0336;%0.028;
RV_Pd_alpha=0.0564;%0.047;


%%%%%%%%%%%%%%%%%%%%%%
%Initialisation of X = PRESSURE FF = FLOW
h=(b-a)/M;
T=zeros(1,M+1);
T=a:h:b;
%Needed for calculation
X(1:15)=0; %X
V(1:15)=0;
V2(1:15)=0;
DX(1:15)=0; %DX/DT
DV(1:15)=0;
FF(1:15)=0; %FLOWS
PP(1:16)=0;
X2(1:15)=0; %X2 (zwischenspeicherung);




%Results Matrices
FLOWS=zeros(15,M+1);
PRESSURES=zeros(15,M+1);
VOLUMES=zeros(15,M+1);
ANGLE=zeros(1,M+1);

%write initial values to array Volumes in ml




X(1)=V0sartUB/C(1);
 X(2)=V0svnUB/C(2);
 X(3)=V0sartLB/C(3); 
 X(4)=V0svnLB/C(4);
 X(7)=V0pa/C(7); 
 X(8)=V0part/C(8); 
 X(9)=V0pvn/C(9); 
 X(12)=V0ao/C(12);
 X(13)=V0lcor/C(13);
 X(14)=0;
 
 V(5)=V0RA;
 V(6)=V0RV;
 V(10)=V0LA;
 V(11)=V0LV;
 
for i=1:15
X2(i)=X(i);
V2(i)=V(i);
end


for j=1:M+1

Q=qheart(T(j));    
ea=eAtrium(T(j));
ev=eVentricle(T(j));

for i=1:15
X(i)=X2(i);
V(i)=V2(i);
end

%Nonlinear Valves
DX(15)=Kav*(PP(11)-X(12))*cos(X(14))-Kp*X(15);
DX(14)=X(15);

%Angle Calculation
 if (X(14)>Angle_max*pi/180) %Saturated Angle Output (Angle max)
 X(14)=Angle_max*pi/180;
 end
  if (X(14)<0)
 X(14)=0;
  end
%Real Crosssection depending on leaflet Angle

A_real=(sin(X(14))/sin(Angle_max/180*pi))^2;


%Equations for CoronaryArteries
FF(13)=(X(13)-0.75*PP(11))/(R(13));
DX(13)=(FF(12)-FF(3)-FF(1)-FF(13))/C(13);

%Equations for SartUB
 FF(1)=(X(1)-X(2))/R(1);
 DX(1)=(FF(12)-FF(3)-FF(13)-FF(1))/(C(1));


%Equations for SvnUB
FF(2)=(X(2)-PP(5))/R(2);
DX(2)=(FF(1)-FF(2))/C(2);

%Equations for SartLB
FF(3)=(X(3)-X(4))/R(3);
DX(3)=(FF(12)-FF(13)-FF(1)-FF(3))/(C(3));

%Equations for SvnLB
FF(4)=(X(4)-PP(5))/R(4);
DX(4)=(FF(3)-FF(4))/C(4);

%Equations for RA
if (PP(5)-PP(6))>0.0001
FF(5)=sqrt(PP(5)-PP(6))*CQtri;
else
FF(5)=0;
end
PP(5)=(V(5))*(ea*(EmaxRA-EdRA)+EdRA);
DV(5)=(FF(4)+FF(2)+FF(13)-FF(5));

% Equations for RV
if (PP(6)-X(7))>0.0001
FF(6)=sqrt(PP(6)-X(7))*CQpa; %No Ventricle FLOW
else
FF(6)=0;
end
PP(6)=(V(6))*ev*EmaxRV+(1-ev)*(RV_Pd_alpha*exp(RV_Pd_kappa*(V(6)))+RV_Pd_beta); %PRESSURE From Elastance
DV(6)=(FF(5)-FF(6));

%Equations for pa
 FF(7)=(X(7)-X(8))/R(7);
 DX(7)=(FF(6)-FF(7))/C(7);

%Equations for part
FF(8)=(X(8)-X(9))/R(8);
DX(8)=(FF(7)-FF(8))/C(8);

%Equations for pvn
FF(9)=(X(9)-PP(10))/R(9);
DX(9)=(FF(8)-FF(9))/C(9);

%Equations for LA
if (PP(10)-PP(11))>0.0001
FF(10)=sqrt(PP(10)-PP(11))*CQmi;
else
FF(10)=0;
end
PP(10)=(V(10))*(ea*(EmaxLA-EdLA)+EdLA);
DV(10)=(FF(9)-FF(10));


%Equations for LV
if (PP(11)-X(12))>0.0001
FF(11)=sqrt(PP(11)-X(12))*CQao*A_real; %No Ventricle FLOW
else
FF(11)=0;
end
PP(11)=(V(11))*ev*EmaxLV+(1-ev)*(LV_Pd_alpha*exp(LV_Pd_kappa*(V(11)))+LV_Pd_beta); %PRESSURE From Elastance
DV(11)=(FF(10)-FF(11)-FF(15));

%Equations for AorticArch
FF(12)=(X(12)-X(3))/(R(12));
DX(12)=(FF(11)-FF(12)+FF(15))/C(12);

%Equations for LVAD
%Inlet Cannula
PP(15)=PP(11)-FF(15)*R(15);
%Outlet Cannula
PP(16)=X(12)+FF(15)*R(16);

%Pump head
%Fitted to pump curve, will be improved when data is available

% FF(15)=(20-(PP(15)-PP(16))*(-0.25))*16.67; 
% 
% if (FF(15)<0)
% FF(15)=0;
% end

FF(15)=0; %LVAD off 


%%%%%%%%%%%%%%%%%%

%Implicit Euler For All compartments except the heart!
X(1)=(X(1)+h*((FF(12)-FF(3)-FF(13))/(C(1))+X(2)/(R(1)*(C(1)))))/(1+h/(R(1)*(C(1)))); %sart UB
X(2)=(X(2)+h*(FF(1)/(C(2))+PP(5)/(R(2)*(C(2)))))/(1+h/(R(2)*(C(2)))); %svnUB
X(3)=(X(3)+h*((FF(12)-FF(13)-FF(1))/(C(3))+X(4)/(R(3)*(C(3)))))/(1+h/(R(3)*(C(3)))); %sartLB
X(4)=(X(4)+h*(FF(3)/C(4)+PP(5)/(R(4)*C(4))))/(1+h/(R(4)*C(4))); %svnLB
X(7)=(X(7)+h*(FF(6)/C(7)+X(8)/(R(7)*C(7))))/(1+h/(R(7)*C(7))); %pa
X(8)=(X(8)+h*(FF(7)/C(8)+X(9)/(R(8)*C(8))))/(1+h/(R(8)*C(8))); %part
X(9)=(X(9)+h*(FF(8)/C(9)+PP(10)/(R(9)*C(9))))/(1+h/(R(9)*C(9))); %pvn
X(12)=(X(12)+h*((FF(11)+FF(15))/C(12)+X(3)/(R(12)*C(12))))/(1+h/(R(12)*C(12))); %AorticArch


%Explicit Euler for Heart
V(5)=V(5)+h*DV(5); %RA
V(6)=V(6)+h*DV(6); %RV

V(10)=V(10)+h*DV(10); %LA
V(11)=V(11)+h*DV(11); %LV


X(13)=X(13)+h*DX(13); %Coronary LEFT

X(14)=X(14)+h*DX(14); %VALVE angle
X(15)=X(15)+h*DX(15); %Derivative of Valve Angle

for i=1:15
X2(i)=X(i);
V2(i)=V(i);
end



for i=1:14  
FLOWS(i,j)=FF(i);
PRESSURES(i,j)=X(i);
VOLUMES(i,j)=X(i)*C(i);
end
ANGLE(1,j)=A_real;

VOLUMES(5,j)=V(5);
PRESSURES(5,j)=PP(5);
VOLUMES(6,j)=V(6);
PRESSURES(6,j)=PP(6);
VOLUMES(10,j)=V(10);
PRESSURES(10,j)=PP(10);
VOLUMES(11,j)=V(11);
PRESSURES(11,j)=PP(11);
FLOWS(11,j)=FF(11);
FLOWS(15,j)=FF(15);
FLOWS(13,j)=(X(12)-0.75*PP(11))/(5*R(13));
end

CO=(mean(FLOWS(11,:))+mean(FLOWS(15,:)))*0.06;%
RES=[VOLUMES(11,:);PRESSURES(11,:)];

plot(RES(1,:),RES(2,:));