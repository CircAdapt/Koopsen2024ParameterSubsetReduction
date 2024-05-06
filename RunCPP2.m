% ---RunCPP2.m---
%
% Run CircAdapt C++ code with limited derivation of model output
% (defined in CircAdaptGetP1.m)

CA = CircAdapt();
CA.buildTriSeg();
CA=setP(CA,P);
CA.runStable();
P=getP1(CA,P);