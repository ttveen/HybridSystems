function [output] = step28fun(input)
x0 = input.x0;
Np = input.Np;
sd_0 = input.sd_0;
sb1_0 = input.sb1_0;
sb2_0 = input.sb2_0;
Ce = input.Ce;
P_load = input.P_load;

%% Load MLD model data
load('Data/MLDmodel.mat')

A = A(1).merged;
B1 = B(1).merged1;
B2 = B(1).merged2;
B3 = B(1).merged3;
B4 = B(1).merged4;

E1 = E(1).merged1;
E2 = E(1).merged2;
E3 = E(1).merged3;
E4 = E(1).merged4;

g5 = g(1).merged;
% x0 = [50; 10; 10];
% P_load = 2*ones(25,1);
%% problem parameters
% Np = 25; % Horizon
Nb = 2; %number of batteries
Wb = [3, 4];
Wd = 10;
Wfuel = 4;
We = 0.4;
T_s = 0.2;

% %
% sd_0 = 0;
% sb1_0 = 0;
% sb2_0 = 0;
%

unknowns = size(B1,2)*Np + size(B2,2)*Np + size(B3,2)*Np;
p_number = 3*Np;
%The V vector are the unknowns, that we optimze for
V_rows = size(B1,2)*Np + size(B2,2)*Np + size(B3,2)*Np + 3*Np + Np; %The 3*Np comes from the p variables, that eliminate the abs in the cost function
M_rows = size(A,1)*Np; %The number of predicted states
P_block = zeros(3*Np,3*Np); %The block that corresponds to the p variables
%% Rewrite the MPC into MILD, co
% Dynamics, the M matrix is the equality constraint on the future states
A_inter = eye(size(A));
M1a = zeros(M_rows,size(B1,2)*Np);
M1b = zeros(M_rows,size(B2,2)*Np);
M1c = zeros(M_rows,size(B3,2)*Np);
M1d = zeros(M_rows,length(B4)); M1d(1:length(B4),1:length(B4)) = eye(length(B4));
for i = 0:Np
    M1a = M1a + kron(diag(ones(1,(25-i)),-i),A_inter*B1);
    M1b = M1b + kron(diag(ones(1,(25-i)),-i),A_inter*B2);
    M1c = M1c + kron(diag(ones(1,(25-i)),-i),A_inter*B3);
    M1d((i+1)*3+1:(i+1)*3+3,:) = A_inter+M1d((i)*3+1:(i)*3+3,:);
    A_inter = A*A_inter;
end
M1d(end-5:end,:) = [];
%Future states do not depend on P, so a zeros block is added, that is
%mulitplied with the p-part of V
M1 = [M1a, M1b, M1c, P_block, zeros(M_rows,Np)];
M2 = [zeros(M_rows-size(A,1),size(A,1));A*A_inter];

% Constraints, F is the linear matrix that we multiply with V, to express
% the linear inequality constraints.
A_inter = eye(size(A));
F1a = zeros(size(E1,1),size(B1,2)*Np); F1a(:,end-2*size(B1,2)+1:end) = [E1*B1, E2];
F1b = zeros(size(E1,1),size(B2,2)*Np); F1b(:,end-2*size(B2,2)+1:end) = [E1*B2, E3];
F1c = zeros(size(E1,1),size(B3,2)*Np); F1c(:,end-2*size(B3,2)+1:end) = [E1*B3, E4];
F1d = zeros(size(E1,1),length(B4));
for i = Np-3:1
    A_inter = A*A_inter;
    F1a(:,i*length(B1,2)+1:3*length(B1,2)) = E1*A_inter*B1;
    F2b(:,i*length(B2,2)+1:3*length(B2,2)) = E1*A_inter*B2;
    F3b(:,i*length(B3,2)+1:3*length(B3,2)) = E1*A_inter*B3;
end
%The inequality constraints on dummy variably p
P_cons = zeros((Np-1)*3*2,V_rows-p_number);
for i = 1:Np
    sd_index(i) = size(B1,2)*Np + i*size(B2,2) - 2;
    sb1_index(i) = size(B1,2)*Np + i*size(B2,2) - 1;
    sb2_index(i) = size(B1,2)*Np + i*size(B2,2) - 0;
    
    p1_index(i) = unknowns + 1 + (i-1)*Np;
    p2_index(i) = unknowns + 2 + (i-1)*Np;
    p3_index(i) = unknowns + 3 + (i-1)*Np;
end
for i = 2:(Np)
    P_cons(2*(i-1)-1,sd_index(i)) = -(-1)^i;
    P_cons(2*(i-1)-1,sd_index(i-1)) = (-1)^i;
    P_cons(2*(i-1),sd_index(i)) = (-1)^i;
    P_cons(2*(i-1),sd_index(i-1)) = -(-1)^i;
end
for i = 2:(Np)
    P_cons(2*(i-1)-1+2*(Np-1),sb1_index(i)) = -(-1)^i;
    P_cons(2*(i-1)-1+2*(Np-1),sb1_index(i-1)) = (-1)^i;
    P_cons(2*(i-1)+2*(Np-1),sb1_index(i)) = (-1)^i;
    P_cons(2*(i-1)+2*(Np-1),sb1_index(i-1)) = -(-1)^i;
end
for i = 2:(Np)
    P_cons(2*(i-1)-1+4*(Np-1),sb2_index(i)) = -(-1)^i;
    P_cons(2*(i-1)-1+4*(Np-1),sb2_index(i-1)) = (-1)^i;
    P_cons(2*(i-1)+4*(Np-1),sb2_index(i)) = (-1)^i;
    P_cons(2*(i-1)+4*(Np-1),sb2_index(i-1)) = -(-1)^i;
end
P_inter =  kron(eye(3),[zeros(2*(Np-1),1), kron(eye((Np-1)),[-1;-1])]);

P_cons = [P_cons,  P_inter];

%constraints on dumy variable p, given some initial condition
P_init_cons = zeros(6,V_rows);
P_init_cons(1,sd_index(1)) = -1; P_init_cons(1,p1_index(1)) = -1;
P_init_cons(2,sd_index(1)) = 1; P_init_cons(1,p1_index(1)) = -1;
P_init_cons(3,sb1_index(1)) = -1; P_init_cons(1,p2_index(1)) = -1;
P_init_cons(4,sb1_index(1)) = 1; P_init_cons(1,p2_index(1)) = -1;
P_init_cons(5,sb2_index(1)) = -1; P_init_cons(1,p3_index(1)) = -1;
P_init_cons(6,sb2_index(1)) = 1; P_init_cons(1,p3_index(1)) = -1;

%Left side of the p constraints
P_cons = [P_cons; P_init_cons];
%RHS
P_rhs = [zeros(2*3*(Np-1),1); ...
    -sd_0; sd_0; -sb1_0; sb1_0; -sb2_0; sb2_0];
%Constraints on P_imp
P_imp_cons = kron(eye(Np),[1,1,1]);
P_imp_cons = [P_imp_cons,zeros(Np,V_rows-size(B1,2)*Np-Np), eye(Np)];

%Merge all the constraints
F1 = [F1a, F1b, F1c, zeros(size(E1,1),p_number), zeros(size(E1,1),Np)];
F1 = [F1;P_cons;P_imp_cons];
F2 = [g5; P_rhs; P_load]; 
F3 = [E1*A_inter; zeros(size(P_cons,1),size(A,1));zeros(size(P_imp_cons,1),size(A,1))];
%% Cost function
% Ce = 50 + 50*sin(pi*T_s*(0:Np-1)/12);
W1 = [zeros(1,unknowns+Np), Wd*ones(1,Np), Wb(1)*ones(1,Np), Wb(2)*ones(1,Np)];
W2 = [Wfuel, We, We, zeros(1,Np*size(A,1)-2*size(A,1)), -Wfuel, -We, -We];
C1 = [zeros(1,V_rows-Np), Ce];
%% Gurobi
%setup

vtype = '';
sense = '';
for i=1:size(B1,2)*Np
   vtype = strcat(vtype,'C');
end
for i=1:size(B2,2)*Np
   vtype = strcat(vtype,'B');
end
for i=1:size(B3,2)*Np
   vtype = strcat(vtype,'C');
end
for i=1:p_number
   vtype = strcat(vtype,'B');
end
for i=1:Np
   vtype = strcat(vtype,'C');
end
for i=1:length(F2)-Np
   sense = strcat(sense,'<'); 
end

for i=1:Np
   sense = strcat(sense,'='); 
end

model.A = sparse(F1);
model.rhs = F2 - F3*x0;
model.vtype = vtype;
model.sense = sense;
model.obj = W1+W2*M1+C1;
% model.lb = zeros(V_rows,1);
model.ub = 200*ones(V_rows,1);
gurobi_write(model, 'mip1.lp');
result = gurobi(model);
output = result.x;


end

