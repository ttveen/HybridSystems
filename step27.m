%% Step 27
clear all
%Load data from step 2.3
load('Data/step23')
%% Constants
Ts = 0.2; %h
eps = 1e-15; %Machine precision for a double
u1 = 5;
u2 = 6.5;
u3 = 11;
%Battery constants
etac = [0.9 0.95];
etad = [0.8 0.77];
xbub = [48 64];
ublb = [-3 -4];
ubub = [2 3];

%Diesel Generator constants
xdlb = 10;
xdub = 100;
udub = 15;
Rf = 0.4;
%% Create the diesel generator
nc = 33; %number of constraints on the generator
nd = 7; %number of deltas
%z inequalities in terms of delta
M = udub;
m = 0;

%inequalities using variable z
zineq = zeros(16,nd);
for i = 1:4
   zineqd(4*(i-1)+1:4*i,i+3) = [-M; m; -m; M];
end
zineqz = zeros(16,4);
for i = 1:4
   zineqz(4*(i-1)+1:4*i,i) = [1; -1; 1; -1];
end
%delta inequalites, result of a product
dproducineq = [blkdiag([0; -1; 1],[0; -1; 1], [0; -1; 1]),...
                blkdiag([1; 1; -1],[1; 1; -1], [1; 1; -1]),...
                repmat([-1;0;1],[3,1])];
%delta constrains, result of bounds on u
duineq = [-blkdiag([u1-M;u1+eps], [u2-M;u2+eps], [u3-M;u3+eps]), zeros(6,4)];
%Create the A matrix
A.d = 1;
%Create the B matrices
B.d1 = 0;
B.d2 = -[zeros(4,1); (par.a1-par.a2); (par.a2-par.a3); (par.a3-par.a4);(par.a4)]';
B.d3 = -[(par.b1 - par.b2); (par.b2 - par.b3); (par.b3 - par.b4); (par.b4)]';
B.d4 = Ts*Rf;
%Create the E and g matrices
E.d1 = [-1; 1; zeros(nc-2,1)];
E.d2 = [zeros(2,1); repmat([0; 0; -1; 1],[4,1]);zeros(9,1); repmat([1;-1],[3,1])];
E.d3 = [zeros(2,nd); zineqd; dproducineq; duineq];
E.d4 = [zeros(2,4); zineqz; zeros(9,4);zeros(6,4)];
g.d = [0;0;repmat([0;0;-m;M],[4,1]);repmat([0;0;1],[3,1]);repmat([M; -eps],[3,1])-[0;u1;0;u2;0;u3]];
%% 2 Batteries
nc = 10; %number of constraints on a battery
%Since the batteries have the same dynamics, but only different parameters,
%a loop can be used to create 2 (but possibly more)
for i = 1:2
    %Create the A matrix
    A(i).b = 1;
    %Create the B matrices
    B(i).b1 = etad(1)*Ts; B(i).b2 = 0; B(i).b3 = (etac(1)-etac(2)); B(i).b4 = 0;
    %Create E and g matices
    E(i).b1 = zeros(nc,1); E(i).b1(1) = -1; E(i).b1(2) = 1;
    E(i).b2 = [0; 0; -1; 1; 1; -1; 0; 0; -1; 1];
    E(i).b3 = [zeros(4,1); -ubub(i); (ublb(i)-eps); -ubub(i); ublb(i); -ublb(i); ubub(i)];
    E(i).b4 = [zeros(6,1); 1; -1; 1; -1];
    g(i).b = [0; xbub(i); -ublb(i); ubub(i); ubub(i); -eps; 0; 0; -ublb(i); ubub(i)];
end
%% Save the data
save('Data/MLDmodel','A','B','E','g');