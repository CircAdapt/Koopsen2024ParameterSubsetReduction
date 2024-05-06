% ---ParSobolObjectiveFunction.m---
%
% Parallel computation of objective function value of Sobol simulations
%
% Koopsen T.
% Last modified: 11/27/2023

parpool(13)

% Patient datasets
datanames = {'MARC1.mat','MARC2.mat','MARC3.mat','MARC4.mat','MARC5.mat','MARC6.mat','MARC7.mat','DEFIMI1.mat','DEFIMI2.mat','DEFIMI3.mat','DEFIMI4.mat','DEFIMI5.mat','DEFIMI6.mat'};

parfor n = 1:13
    nsuc = SobolObjectiveFunction(datanames{n});
end