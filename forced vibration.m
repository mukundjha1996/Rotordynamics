%assignMent 4
%Date of birth Month is NoveMber

%question no. 31 is not available in exercise hence question  no. 12 is
%solveD

clc
clear all
syMs l1 l2 l3 l4 l5 g w D1 D2 D3 D4 D

J=pi*D^4/32;
Disp('............FINITE ELEMENT METHOD.........')

M=cell(1,5);
K=cell(1,5);

%eleMent 1
IP1=(4*0.04^2)/2;
M{1}=[0 0;0 IP1];
K1=(g*J/l1)*[1 -1;-1 1];
K{1}=Double(subs(K1,{g,D,l1},{0.8*10^11,0.02,0.15}));

%eleMent 2
IP2=(5*0.05^2)/2;
M{2}=[0 0;0 IP2];
K2=(g*J/l2)*[1 -1;-1 1];
K{2}=Double(subs(K2,{g,D,l2},{0.8*10^11,0.02,0.05}));

%eleMent 3
IP3=(6*0.06^2)/2;
M{3}=[0 0;0 IP3];
K3=(g*J/l3)*[1 -1;-1 1];
K{3}=Double(subs(K3,{g,D,l3},{0.8*10^11,0.02,0.05}));

%eleMent 4
IP4=(7*0.07^2)/2;
M{4}=[0 0;0 IP4];
K4=(g*J/l4)*[1 -1;-1 1];
K{4}=Double(subs(K4,{g,D,l4},{0.8*10^11,0.02,0.05}));

%eleMent 5
M{5}=[0 0;0 0];
K5=(g*J/l5)*[1 -1;-1 1];
K{5}=Double(subs(K5,{g,D,l5},{0.8*10^11,0.02,0.15}));


%global Mass Matrix
MASS_MATRIX=zeros(6,6);
for i=1:5
    MASS_MATRIX(i:i+1,i:i+1)=MASS_MATRIX(i:i+1,i:i+1)+M{i};
enD
  
%global stiffness Matrix
STIFFNESS_MATRIX=zeros(6,6);
for i=1:5
    STIFFNESS_MATRIX(i:i+1,i:i+1)=STIFFNESS_MATRIX(i:i+1,i:i+1)+K{i};
enD
MASS_MATRIX
STIFFNESS_MATRIX
%bounDary conDition
MASS_MATRIX=vpa(-w^2*MASS_MATRIX);
A2=vpa(MASS_MATRIX+STIFFNESS_MATRIX);


%solving for frequency
A2=Det(A2);
frequency2=solve(A2,w);
Disp('THE FREQUENCY USING FINITE ELEMENT METHOD IS:')
frequency=frequency2(5:8)
Disp('------------------------------------------')
