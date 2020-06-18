x = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)';
a = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
y = a*x+(-4+8*rand(length(a),1));

x_hat = intvar(6,1);

residuals = y-a*x_hat;
bound = sdpvar(length(residuals),1);
F = [-bound <= residuals <= bound];
optimize(F,sum(bound));
x_L1 = value(x_hat);
optimize([],residuals'*residuals);
x_L2 = value(x_hat);
bound = sdpvar(1,1);
F = [-bound <= residuals <= bound];
optimize(F,bound);
x_Linf = value(x_hat);