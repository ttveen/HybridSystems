%% Step 29
clear all
load('Data/MLDmodel.mat')

A = A(1).merged;
B1 = B(1).merged1;
B2 = B(1).merged2;
B3 = B(1).merged3;
B4 = B(1).merged4;
%

T_s = 0.2;
Np= 25;
input.sd_0 = 0;
input.sb1_0 = 0;
input.sb2_0 = 0;

input.Ce = 50 + 50*sin(pi*T_s*(0:Np-1)/12);
input.x0 = [50;10;10];
input.Np = 25;

Tfinal = 100;
P_load = zeros(Tfinal+Np,1); 
P_load(21:50) = 30+2*(21:50);
P_load(51:Tfinal+Np) = 45;
input.P_load = P_load(1:Np);
x = input.x0;
for k=1:Tfinal
    input.x0 = x(:,k)
    output = step28fun(input);
    input.sd_0 = output(3*Np + 9 - 2);
    input.sb1_0 = output(3*Np + 9 - 1);
    input.sb2_0 = output(3*Np + 9 - 0);
    input.P_load = P_load(k:k+Np-1); 
    input.Ce = 50 + 50*sin(pi*T_s*(k:k+Np-1)/12);
%     input.x0 = output(1:3);
    
    x(:,k+1) = A*x(:,k)+B1*output(1:3)+B2*output(3*Np+1:3*Np+9)+B3*output(3*Np+9*Np+1:3*Np+9*Np+6)+B4;
    ud(k) = output(1);
    ub1(k) = output(2);
    ub2(k) = output(2);
    P_imp(k) = output(end-24);
end