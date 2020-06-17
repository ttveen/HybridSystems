%% Setup
clear; clc; close all;

%% Load MLD model data
load('Data/MLDmodel.mat')

%% MILP Yalmip

% problem parameters
Np = 25; % Horizon
Nb = 2; %number of batteries
Wb = [3, 4];
Wd = 10;
Wfuel = 4;
We = 0.4;

% switching signals
sb = binvar(repmat(2,1,Np),repmat(1,1,Np));
delta = binvar(repmat(7,1,Np),repmat(1,1,Np));

sb_0 = binvar(repmat(2,1,1),repmat(1,1,1));
sd_0 = binvar(repmat(1,1,1),repmat(1,1,1));

% dummy variable
p = binvar(repmat(3,1,Np),repmat(1,1,Np));

% sdpvar(repmat(nu,1,N),repmat(1,1,N));
% batteries 
xb = sdpvar(repmat(2,1,Np+1),repmat(1,1,Np+1));
ub = sdpvar(repmat(2,1,Np),repmat(1,1,Np));
zb = sdpvar(repmat(2,1,Np),repmat(1,1,Np));

% diesel generators
xd = sdpvar(repmat(1,1,Np+1),repmat(1,1,Np+1));
ud = sdpvar(repmat(1,1,Np),repmat(1,1,Np));
zd = sdpvar(repmat(4,1,Np),repmat(1,1,Np));

% imported power
P_imp = sdpvar(repmat(1,1,Np),repmat(1,1,Np));
P_load = sdpvar(repmat(1,1,Np),repmat(1,1,Np));
Ce = sdpvar(repmat(1,1,Np),repmat(1,1,Np));



objective = 0;
constraints = [];

for j = 1:Np
    
    objective = objective + Wb(1)*p{j}(1) + Wb(2)*p{j}(2) + Wd*p{j}(3) ...
        + P_imp{j}*Ce{j};
    
    % diesel generator dynamics
    constraints = [constraints; xd{j+1} == A(1).d*xd{j}+B(1).d1*ud{j}...
        + B(1).d2*delta{j} + B(1).d3*zd{j} + B(1).d4];
    
    % diesel generator inequalities
    constraints = [constraints; E(1).d1*xd{j}+E(1).d2*ud{j}...
        + E(1).d3*delta{j} + E(1).d4*zd{j} <= g(1).d];
    
    % batteries dynamic
    % battery 1
    constraints = [constraints; xb{j+1}(1) == A(1).b*xb{j}(1)+B(1).b1*ub{j}(1)...
        + B(1).b2*sb{j}(1) + B(1).b3*zb{j}(1) + B(1).b4];
    
    % battery 2
    constraints = [constraints; xb{j+1}(2) == A(2).b*xb{j}(2)+B(2).b1*ub{j}(2)...
        + B(2).b2*sb{j}(2) + B(2).b3*zb{j}(2) + B(2).b4];
    
    
    % batteries inequalities
    % battery 1
    constraints = [constraints; E(1).b1*xb{j}(1) + E(1).b2*ub{j}(1)...
        + E(1).b3*sb{j}(1) + E(1).b4*zb{j}(1) <= g(1).b];
    
    % battery 2
    constraints = [constraints; E(2).b1*xb{j}(2) + E(2).b2*ub{j}(2)...
        + E(2).b3*sb{j}(2) + E(2).b4*zb{j}(2) <= g(2).b];
   
    
end

for j = 2:Np
    constraints = [constraints; -p{j}(1:2) <= sb{j}-sb{j-1} <= p{j}(1:2)];
    constraints = [constraints; -p{j}(3) <= delta{j}(end)-delta{j-1}(end) <= p{j}(3)];
end

constraints = [constraints; -p{j}(1:2) <= sb{j}-sb_0 <= p{j}(1:2)];
constraints = [constraints; -p{j}(3) <= delta{j}(end)-sd_0 <= p{j}(3)];

objective = objective - Wfuel*(xd{Np+1}-xd{1}) -We*sum(xb{Np+1}-xb{1});

%% create controller

parameters_in = {xd{1},xb{1},sd_0,sb_0,[P_load{:}],[Ce{:}]};

solutions_out = {[xd{:}],[xb{:}],[ud{:}],[ub{:}],[delta{:}],[sb{:}],[P_imp{:}]};

%The controller, given the constraints and objective
options = sdpsettings('verbose',1,'solver','gurobi');
controller = optimize(constraints,objective,options,parameters_in,solutions_out)

%% Simulate
Tfinal = 200;
T_s = 0.20;
% t = 0:T_s:Tfinal;

xd_ = 50; 
xb_ = [10; 10];
sd_0_ = 0;
sb_0_ = [0; 0];

P_load_ = zeros(1,Tfinal+Np); 
P_load_(21:50) = 30+2*(21:50);
P_load_(51:Tfinal+Np) = 45;

for k = 1:Tfinal
    Ce_ = 50+50*sin((pi*T_s*(k:k+Np-1))/12);
    
    output = controller(xd_,xb_,sd_0_,sb_0_,P_load_(k:k+Np-1),Ce_); 
end
    
    
