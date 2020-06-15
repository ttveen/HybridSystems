%% Step 2.3
% The objective is to fit 4 linear functions, parameterised with a_i and
% b_i, to a piecewise function, consisting of 5 nonlinear functions. The
% cost function is given in the exercise. It is the integral of the
% nonlinear function minus the linear fit squared.
%
% We have minimised the cost as follows:
% First approximate the first part, results in a_1 and b_1.
% Then, compute f(u_d) at the end of the first part (at u_d = 5)
% Solve for the second part (5<u_d<6.5) with an equality constraint a1+b1*5
% = a2+b2*5.
% Repeat this for all parts.
%
% At the end of the file, all functions are created. The represent the cost
% for each piecewise part, partioned like the affine fit.
%
% For all minimisations, fmcon() is used
%
% A plot of the nonlinear and affine piecewise function is added to
% compare.
close all
clear all

%Values that determine the parts of the affine function
u1 = 5;
u2 = 6.5;
u3 = 11;

%Find the first approximation, a1 and b1
x0 = [1,1]; %initial condition
f1hat = fmincon(@funf1,x0);
par.a1 = f1hat(1,1); par.b1 = f1hat(1,2);
%Compute the boundary at u_d = 5 = u1
Aeq = [1 u1];
beq = par.a1+par.b1*u1;

%The fit to second part should start at the right boundary condition
f2hat = fmincon(@funf2,[-80,20],[],[],Aeq,beq);
par.a2 = f2hat(1,1); par.b2 = f2hat(1,2);

%The third part also should begin at the right boundary
Aeq = [1 u2];
beq = par.a2+par.b2*u2;
f3hat = fmincon(@funf3,x0,[],[],Aeq,beq);
par.a3 = f3hat(1,1); par.b3 = f3hat(1,2);

%Same for fourth section
Aeq = [1 u3];
beq = par.a3+par.b3*u3;
f4hat = fmincon(@funf4,[-200,50],[],[],Aeq,beq);
par.a4 = f4hat(1,1); par.b4 = f4hat(1,2);



%% Create plots
%First create the nonlinear function
syms ud
% assume(0<= ud <=15)
f1 = piecewise(0<=ud<2,ud^2+4,...
               2<=ud<5, 4*ud,...
               5<=ud<7, -9.44*ud^3+166.06*ud^2-948.22*ud+1790.28,...
               7<=ud<9, -11.78*ud+132.44,...
               9<=ud<=15, 4.01*(ud-10.47)^2+17.79);

%Create the affine function
syms a1 a2 a3 a4
syms b1 b2 b3 b4

u1 = 5;
u2 = 6.5;
u3 = 11;
f2(a1,a2,a3,a4,b1,b2,b3,b4) = piecewise(0<=ud<u1, a1+b1*ud,...
                    u1<=ud<u2, a2+b2*ud,...
                    u2<=ud<u3, a3+b3*ud,...
                    u3<=ud<=15, a4+b4*ud);


%% plot
figure
hold on
fplot(f1,[0 15])
fplot(f2(par.a1,par.a2,par.a3,par.a4,par.b1,par.b2,par.b3,par.b4),[0 15])
hold off
xlabel({'$u_d(k)$'},'Interpreter', 'latex')
ylabel({'function value'},'Interpreter', 'latex')
legend({'$f(u_d(k))$','$\hat{f}(u_d(k))$'},'Interpreter', 'latex')
step23plot = gcf;
%Save the plot
saveas(step23plot,'Latex/images/step23','eps')
%Save the parameters
save('Data/step23.mat', 'par')
%% Functions
% There is a function that returns the cost as function of x = [a_i,b_i].
% fi reppresents the cost of the part i, of the affine piecewise function.
function f1 = funf1(x)
    a1 = x(1); b1 = x(2);
    fun1 = @(u,a1,b1) (u.^2 + 4 - (a1+b1*u)).^2;
    fun2 = @(u,a1,b1) (4*u - (a1+b1*u)).^2;
    f1 = integral(@(u) fun1(u,a1,b1),0,2) + integral(@(u) fun2(u,a1,b1),2,5);
end
function f2 = funf2(x)
    a2 = x(1); b2 = x(2);
    fun1 = @(u,a2,b2) ((-9.44*u.^3 + 166.06*u.^2 - 948.22*u + 1790.28) - (a2+b2*u)).^2;
    f2 = integral(@(u) fun1(u,a2,b2),5,6.5);
end
function f3 = funf3(x) 
    a3 = x(1); b3 = x(2);
    fun1 = @(u,a3,b3) ((-9.44*u.^3 + 166.06*u.^2 - 948.22*u + 1790.28) - (a3+b3*u)).^2;
    fun2 = @(u,a3,b3) ((-11.78*u + 132.44) - (a3+b3*u)).^2;
    fun3 = @(u,a3,b3) ((4.01*(u -10.47) + 17.79) - (a3+b3*u)).^2;
    f3 = integral(@(u) fun1(u,a3,b3),6.5,7) + integral(@(u) fun2(u,a3,b3),7,9) + integral(@(u) fun3(u,a3,b3),9,11);
end
function f4 = funf4(x)
    a4 = x(1); b4 = x(2);
    fun1 = @(u,a4,b4) ((4.01*(u -10.47).^2 + 17.79) - (a4+b4*u)).^2;
    f4 = integral(@(u) fun1(u,a4,b4),11,15);
end