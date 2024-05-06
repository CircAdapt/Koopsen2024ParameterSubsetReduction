% ---RunCPP1.m---
%
% Run CircAdapt C++ code with limited definition of state variables 
% (defined in CircAdaptSetP1.m) and limited derivation of model output 
% (defined in CircAdaptGetP1.m)

CA=setP1(CA,P);
CA.runStable();
P=getP1(CA,P);