% ---RunCPP.m---
%
% Run CircAdapt C++ code

CA = CircAdapt();
CA.buildTriSeg();

CA=setP(CA,P);
CA.setScalar('','','','solverTol',1e-4)
CA.setScalar('','','','minTol',1e-3)
CA.setScalar('','','','maxBeats',100)

if P.General.PressFlowContr == 1
    CA.runStable();
else
    for beat = 1:5
        CA.runSingleBeat();
    end
end
P=getP(CA);