% ---SobolObjectiveFunction.m---
%
% Derive matrix with objective function values of Sobol simulations
%
% Koopsen T.
% Last modified: 11/27/2023

function nsuc = SobolObjectiveFunction(data)

% Input arguments:
% data    - data struct containing patient LV volume and strain data
%
% Output arguments:
% nsuc    - number of successful simulations (computationally stable)

nsubdata = 30; % number of subdata matrices
nindex = 100000; % number of simulations per matrix

% Objective function has 15 components: Chi2VED, Chi2EF, Chi2tcycle, Chi2Ei(*6 segments), Chi2ERi (*6 segments)
all_OF = zeros(nindex,nsubdata,15)+Inf;

load(data,'measurement')

nsuc = 0; % number of successful simulations (computationally stable)

% Fill objective function value matrix
for n = 1:nsubdata
    load(['Output_Diaphony_sub',num2str(n),'.mat'],['output_model',num2str(n)])
    eval(['dataname = output_model',num2str(n),';'])
    for i = 1:nindex
        if isstruct(dataname{i})
            nsuc = nsuc+1;
            D = dataname{i};
            OF_val = CalculateObjectiveFunction(D,measurement);
            for k = 1:15
                all_OF(i,n,k) = OF_val(k);
            end
        end
    end
    clear(['output_model',num2str(n)])
end

save(['OF_',data],'all_OF')
