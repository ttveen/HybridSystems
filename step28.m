%% Setup
clear; close all;

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
T_s = 0.2;

yalmip('clear')
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
% P_load = 50;
% Ce = 10;
% k = sdpvar;

% objective = xb{1}% + ud{1} + delta{1}(1);
% constraints = [-3 <= xb{1} <= 3];
constraints = [];
objective = 0;
objective = objective - Wfuel*(xd{Np+1}-xd{1}) -We*sum(xb{Np+1}-xb{1});

for j = 1:Np
%     Ce = 50+50*sin((pi*T_s*(k))/12);
    
    objective = objective + Wb(1)*p{j}(1) + Wb(2)*p{j}(2) + Wd*p{j}(3) ...
        + P_imp{j}*Ce{j};
    
    % diesel generator dynamics
    constraints = [constraints; xd{j+1} == A(1).d*xd{j}+B(1).d1*ud{j}...
        + B(1).d2*delta{j} + B(1).d3*zd{j} + B(1).d4];
    
    % diesel generator inequalities !!!!!!!!!!!!!
    constraints = [constraints; E(1).d1*xd{j}+E(1).d2*ud{j}...
        + E(1).d3*delta{j} + E(1).d4*zd{j} <= g(1).d];
%     
%     % batteries dynamic
%     % battery 1
    constraints = [constraints; xb{j+1}(1) == A(1).b*xb{j}(1)+B(1).b1*ub{j}(1)...
        + B(1).b2*sb{j}(1) + B(1).b3*zb{j}(1) + B(1).b4];
%     
%     % battery 2
    constraints = [constraints; xb{j+1}(2) == A(2).b*xb{j}(2)+B(2).b1*ub{j}(2)...
        + B(2).b2*sb{j}(2) + B(2).b3*zb{j}(2) + B(2).b4];
%     
%     
%     % batteries inequalities
%     % battery 1
    constraints = [constraints; E(1).b1*xb{j}(1) + E(1).b2*ub{j}(1)...
        + E(1).b3*sb{j}(1) + E(1).b4*zb{j}(1) <= g(1).b];
%     
%     % battery 2
    constraints = [constraints; E(2).b1*xb{j}(2) + E(2).b2*ub{j}(2)...
        + E(2).b3*sb{j}(2) + E(2).b4*zb{j}(2) <= g(2).b];
%    
    constraints = [constraints; P_imp{j} == P_load{j} - ud{j} - sum(ub{j})];
end

for j = 2:Np
    constraints = [constraints; -p{j}(1:2) <= sb{j}-sb{j-1} <= p{j}(1:2)];
    constraints = [constraints; -p{j}(3) <= delta{j}(end)-delta{j-1}(end) <= p{j}(3)];
end

constraints = [constraints; -p{1}(1:2) <= sb{1}-sb_0 <= p{1}(1:2)];
constraints = [constraints; -p{1}(3) <= delta{1}(end)-sd_0 <= p{1}(3)];



%% create controller

parameters_in = {xd{1},xb{1},sd_0,sb_0,[P_load{:}],[Ce{:}]};

solutions_out = {[xd{:}],[xb{:}],[ud{:}],[ub{:}],[delta{:}],[sb{:}],[P_imp{:}]};

%The controller, given the constraints and objective
options = sdpsettings('verbose',1,'solver', 'gurobi','debug',1);
controller = optimizer(constraints,objective,options,parameters_in,solutions_out);
%%
%% Simulate
Tfinal = 200;
T_s = 0.20;
% t = 0:T_s:Tfinal;

xd = zeros(1,Tfinal+1);
xb = zeros(2,Tfinal+1);
ud = zeros(1,Tfinal);
ub = zeros(2,Tfinal);
P_imp = zeros(1,Tfinal);


xd(1) = 50; 
xb(:,1) = [10; 10];
sd_0 = 0;
sb_0 = [0; 0];

P_load = zeros(1,Tfinal+Np); 
P_load(21:50) = 30+2*(21:50);
P_load(51:Tfinal+Np) = 45;

for k = 1:Tfinal
%     Ce = 50+50*sin((pi*T_s*(k-1:k+Np-2))/12);
    Ce = zeros(1,25);
    inputs = {xd(k),xb(:,k),sd_0,sb_0,P_load(k:k+Np-1),Ce};
    output = controller{inputs}; 
    xd(k+1) = output{1}(2);
    xb(:,k+1) = output{2}(:,2);
    ud(k+1) = output{3}(1);
    ub(:,k+1) = output{4}(:,1);
    sd_0 = output{5}(end,1);
    sb_0 = output{6}(:,1);
    P_imp(k) = output{7}(1);
    
    disp(k)
end
    
plot(ud)
