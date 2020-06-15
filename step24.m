%% Step 2.4
close all
clear all
load('Data/step23.mat')
x0 = [par.a1 par.a2 par.a3 par.a4 par.b1 par.b2 par.b3 par.b4 5 6.5 11];
parmin = fmincon(@funcost,x0,[],[],[],[],[],[],@nonlincon);
par24.a1 = parmin(1); par24.a2 = parmin(2); par24.a3 = parmin(3); par24.a4 = parmin(4);
par24.b1 = parmin(5); par24.b2 = parmin(6); par24.b3 = parmin(7); par24.b4 = parmin(8);
par24.u1 = parmin(9); par24.u2 = parmin(10); par24.u3 = parmin(11);
%% Plot
syms ud
f1 = piecewise(0<=ud<2, ud^2+4,...
           2<=ud<5, 4*ud,...
           5<=ud<7, -9.44*ud^3+166.06*ud^2-948.22*ud+1790.28,...
           7<=ud<9, -11.78*ud+132.44,...
           9<=ud<=15, 4.01*(ud-10.47)^2+17.79);
f2 = piecewise(0<=ud<par24.u1, par24.a1+par24.b1*ud,...
   par24.u1<=ud<par24.u2, par24.a2+par24.b2*ud,...
   par24.u2<=ud<par24.u3, par24.a3+par24.b3*ud,...
   par24.u3<=ud<=15, par24.a4+par24.b4*ud);
figure()
hold on
fplot(f1,[0 15])
fplot(f2,[0 15])
hold off
xlabel({'$u_d(k)$'},'Interpreter', 'latex')
ylabel({'function value'},'Interpreter', 'latex')
legend({'$f(u_d(k))$','$\hat{f}(u_d(k))$'},'Interpreter', 'latex')
%% Save the data
step24plot = gcf;
%Save the plot
saveaspdf(step24plot,'Latex/images/step24')
%Save the parameters
save('Data/step24.mat', 'par')
%% Functions
function cost = funcost(x)
    syms ud
    f1 = piecewise(0<=ud<2, ud^2+4,...
               2<=ud<5, 4*ud,...
               5<=ud<7, -9.44*ud^3+166.06*ud^2-948.22*ud+1790.28,...
               7<=ud<9, -11.78*ud+132.44,...
               9<=ud<=15, 4.01*(ud-10.47)^2+17.79);
           
    a1 = x(1); a2 = x(2); a3 = x(3); a4 = x(4);
    b1 = x(5); b2 = x(6); b3 = x(7); b4 = x(8);
    u1 = x(9); u2 = x(10); u3 = x(11);
   f2 = piecewise(0<=ud<u1, a1+b1*ud,...
       u1<=ud<u2, a2+b2*ud,...
       u2<=ud<u3, a3+b3*ud,...
       u3<=ud<=15, a4+b4*ud);
   
   cost = double(int((f1-f2)^2,ud,0,15));
end
function [c,ceq] = nonlincon(x)
    c = [];
    ceq = [x(1)+x(5)*x(9) - (x(2)+x(6)*x(9));
           x(2)+x(6)*x(10) - (x(3)+x(7)*x(10));
           x(3)+x(7)*x(11) - (x(4)+x(8)*x(11))];
end