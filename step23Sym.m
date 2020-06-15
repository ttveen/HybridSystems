%% Setup
clear; clc; close all;


%% Create the piecewise functions

% non-linear PWA function of the diesel generator's fuel consumption
syms ud

f1 = piecewise(0<=ud<2, ud^2+4,...
               2<=ud<5, 4*ud,...
               5<=ud<7, -9.44*ud^3+166.06*ud^2-948.22*ud+1790.28,...
               7<=ud<9, -11.78*ud+132.44,...
               9<=ud<=15, 4.01*(ud-10.47)^2+17.79);

% PWA approximation of the above PWA function
syms a [1 4]
syms b [1 4]
u1 = 5;
u2 = 6.5;
u3 = 11;

f2(a,b) = piecewise(0<=ud<u1, a1+b1*ud,...
                    u1<=ud<u2, a2+b2*ud,...
                    u2<=ud<u3, a3+b3*ud,...
                    u3<=ud<=15, a4+b4*ud);
               
                
%% Determine a and b values

% create objective function                
f = int((f1-f2)^2,ud,0,15);

% calculate the partial derivatives of f over a and b variables 
diff_eqn = jacobian(f,[a b]);

% set the partial derrivatives equal to zero and solve for a and b 
sol = solve(diff_eqn == 0,[a b]);

% assign solution to individual variables
s_a1 = sol.a1;
s_a2 = sol.a2;
s_a3 = sol.a3;
s_a4 = sol.a4;

s_a = double([s_a1; s_a2; s_a3; s_a4]);

s_b1 = sol.b1;
s_b2 = sol.b2;
s_b3 = sol.b3;
s_b4 = sol.b4;

s_b = double([s_b1; s_b2; s_b3; s_b4]);
%% Plot the AWP comparison
figure ('Name','AWP Comparison')
xlabel('xlabel');
ylabel('ylabel');
hold on
fplot(f1,[0 15])
fplot(f2(s_a1,s_a2,s_a3,s_a4,s_b1,s_b2,s_b3,s_b4),[0 15])
hold off

saveaspdf(gcf,'Latex/images/part23Sym.pdf')