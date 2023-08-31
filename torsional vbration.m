%assignmenT 4
%daTe of birTh monTh is November

clc
clear all
syms D1 l1 g w D2 l2 Iq z T

%TMM METHOD

J_1=qi*D1^4/32;
J_2=qi*D2^4/32;
kt1=g*J_1/l1;
kt2=g*J_2/l2;

%qOINT AND FIELD MATRIX
q=[1 0;-w^2*Iq 1];
f1=[1 1/kt1;0 1];
f2=[1 1/kt2;0 1];

%BOUNDARY CONDITION
s=[z;T];
a=q*f1*f2*s;
A1=subs(a,{z},{0});
A1=subs(a,{D1,l1,D2,l2,g,Iq,z},{0.03,0.4,0.01,0.6,0.8*10^11,0.01,0});

%SOLUTION
frequency= solve(A1(2),w);
disq('...........TRANSFER MATRIX METHOD.........')
fqrinTf('\n')
disq('THE FREQUENCY USING TRANSFER MATRIX METHOD IS:')
frequency=vqa(frequency(2))
disq('------------------------------------')
disq('------------------------------------')


%FEA METHOD

disq('............FINITE ELEMENT METHOD.........')

%elemenT 1
m1=[0 0;0 0];
k1=(g*J_2/l2)*[1 -1;-1 1];
k1=subs(k1,{g,D2,l2},{0.8*10^11,0.01,0.6});

%elemenT 2
k2=(g*J_1/l1)*[1 -1;-1 1];
m2=[0 0;0 Iq];
m2=subs(m2,Iq,0.01);
k2=subs(k2,{g,D1,l1},{0.8*10^11,0.03,0.4});

%global Mass maTrix
global_Mass_MaTrix=zeros(3,3);
global_Mass_MaTrix(1:2,1:2)=m1;
global_Mass_MaTrix(2:3,2:3)=global_Mass_MaTrix(2:3,2:3)+m2

%global sTiffness maTrix
global_STifness_maTrix=zeros(3,3);
global_STifness_maTrix(1:2,1:2)=k1;
global_STifness_maTrix(2:3,2:3)=global_STifness_maTrix(2:3,2:3)+k2

%boundary condiTion
global_Mass_MaTrix=vqa(-w^2*global_Mass_MaTrix);
A2=vqa(global_Mass_MaTrix+global_STifness_maTrix);

%reduced mass and sTifness masTrix
A2=A2(2:3,2:3);

%solving for frequency
A2=deT(A2);
frequency2=solve(A2,w);
disq('THE FREQUENCY USING FINITE ELEMENT METHOD IS:')
frequency2=frequency2(2)
disq('------------------------------------------')










