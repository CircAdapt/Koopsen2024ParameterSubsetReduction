% ---CalculateDiaphony.m---
%
% Calculate parameter diaphony
%
% Koopsen T.
% Last modified: 11/28/2023

function diaphony = CalculateDiaphony(index,input)

% Input arguments:
% index    - indices of best simulations in input matrix
% input    - input matrix of Sobol simulations

% Matrix with parameters of best N simulations
all_pars = input(index,:);

% Calculate diaphony for all parameters
sigmaN = zeros(103,1);

for n = 1:length(index)
    sigmaN = sigmaN + exp(1i.*2.*pi.*all_pars(n,:)');
end
sigmaNnorm = sigmaN./length(index);
diaphony = abs(sigmaNnorm);