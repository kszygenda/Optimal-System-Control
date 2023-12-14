clc; close all;
% Parametry i macierze
syms s k1 k2 k3
K=[k1 k2 k3];
c1=1;
c2=2;
c3=3;
A=[-1/c1 0 0; 0 -1/c2 0;0 0 -1/c3];
B=[1/c1;1/c2;1/c3];
C=[1 0 0];
D=0;
sys=ss(A,B,C,D);

% Tutaj macierz sterowalna
P=[B A*B A*A*B];
As=P*A*P^(-1);
Bs=P*B;
Cs=C*P^(-1);
rownanie=det(s*eye(3,3) - (As-Bs*K))==(s+1)*(s+2)*(s+5);
solve(rownanie,[k1 k2 k3])
rownania=[5/6*k1+5/12*k2+47/162*k3+11/6==8, 2/3*k1+5/16*k2+2/9*k3+1==17,1/6*k1+1/16*k2+7/162*k3+1/6==10];
solve(rownania,[k1 k2 k3])

