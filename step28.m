%% Step 28
load('Data/MLDmodel')
%Constants
Nb = 2; %number of batteries
Wb = [3, 4];
Wd = 10;
Wfuel = 4;
We = 0.4;
%%
names = {'x';'y';'z'};
model.A = sparse([1 2 3; 1 1 0]);
model.obj = [1 1 2];
model.rhs = [4; 1];
model.sense = '<>';
model.vtype = 'B';
model.modelsense = 'max';
model.varnames = names;
gurobi_write(model, 'mip1.lp');

params.outputflag = 0;

result = gurobi(model, params);

disp(result);

for v=1:length(names)
    fprintf('%s %d\n', names{v}, result.x(v));
end

fprintf('Obj: %e\n', result.objval);