% ---setSim.m---
%
% Set parameter values of simulation
%
% Koopsen T.
% Last modified: 11/27/2023

function P = setSim(refP,par_set,pars)

% Input arguments:
% refP    - reference P-struct
% par_set - parameter set used for optimization
% pars    - array(1*270) with normalized parameter values, indexed as follows:
% 1        - q0
% 2        - p0
% 3        - k_Sy
% 4        - p0_Sy
% 5        - A0_Sy
% 6        - dTauAv
% 7-24     - dT_Lv (18 patches)
% 25       - dT_La
% 26-43    - VWall_Lv (18 patches)
% 44       - VWall_La
% 45-62    - AmRef_Lv (18 patches)
% 63       - AmRef_Rv
% 64       - AmRef_La
% 65-82    - SfAct_Lv (18 patches)
% 83       - SfAct_Rv
% 84       - SfAct_La
% 85-102   - vMax_Lv (18 patches)
% 103-120  - TR_Lv (18 patches)
% 121      - TR_Rv
% 122      - TR_La
% 123-140  - TD_Lv (18 patches)
% 141-158  - Ls0Pas_Lv (18 patches)
% 159      - Ls0Pas_La
% 160-177  - SfPas_Lv (18 patches)
% 178-195  - k1_Lv (18 patches)
% 196      - k1_La
% 197-214  - dLsPas_Lv (18 patches)
% 215      - LenSE_Lv_G (global)
% 216-233  - ADO_Lv (18 patches)
% 234      - ADO_La
% 235-252  - LDAD_Lv (18 patches)
% 253-270  - LDCI_Lv (18 patches)

% Load reference P-struct
load(refP,'P');

% Global circulation parameters (3)
P.General.q0 = pars(1);
P.General.p0 = pars(2);
%P.General.tCycle = pars(271);

% ArtVen parameters (3)
P.ArtVen.k(:,1) = [8;10].*pars(3);
P.ArtVen.p0(:,1) = [1.2166e4;0.0141e4].*pars(4);
P.ArtVen.A0(:,1) = [0.4971e-3;0.4992e-3].*pars(5);

% Electrical substrate (20)
P.General.dTauAv = pars(6);
if contains('dT_Lv',par_set)
    P.Patch.dT(3:20) = pars(7:24);
elseif contains('dT_g_Lv',par_set)
    P.Patch.dT(3:20) = repmat(pars(7),1,18);
end
P.Patch.dT(1) = pars(25);

% Other patch parameters (245)
% VWall (19)
if contains('VWall_Lv',par_set)
    P.Patch.VWall(3:20) = pars(26:43);
end
if contains('VWall_g_Lv',par_set)
    P.Patch.VWall(3:14) = repmat(pars(26)*8.0007e-6,1,12);
    P.Patch.VWall(15:20) = repmat(pars(26)*5.3926e-6,1,6);
end
P.Patch.VWall(1) = pars(44);
% AmRef (20)
if contains('AmRef_Lv',par_set)
    P.Patch.AmRef(3:20) = pars(45:62);
end
if contains('AmRef_g_Lv',par_set)
    P.Patch.AmRef(3:14) = repmat(pars(45)*8.1711e-4,1,12);
    P.Patch.AmRef(15:20) = repmat(pars(45)*8.1404e-4,1,6);
end
P.Patch.AmRef(21) = pars(63);
P.Patch.AmRef(1) = pars(64);
% SfAct (20)
if contains('SfAct_Lv',par_set)
    P.Patch.SfAct(3:20) = pars(65:82);
end
if contains('SfAct_g_Lv',par_set)
    P.Patch.SfAct(3:20) = repmat(pars(65),1,18);
end
P.Patch.SfAct(21) = pars(83);
P.Patch.SfAct(1) = pars(84);
% vMax (18)
if contains('vMax_Lv',par_set)
    P.Patch.vMax(3:20) = pars(85:102);
end
if contains('vMax_g_Lv',par_set)
    P.Patch.vMax(3:20) = repmat(pars(85),1,18);
end
% TR (20)
if contains('TR_Lv',par_set)
    P.Patch.TR(3:20) = pars(103:120);
end
if contains('TR_g_Lv',par_set)
    P.Patch.TR(3:20) = repmat(pars(103),1,18);
end
P.Patch.TR(21) = pars(121);
P.Patch.TR(1) = pars(122);
% TD (18)
if contains('TD_Lv',par_set)
    P.Patch.TD(3:20) = pars(123:140);
end
if contains('TD_g_Lv',par_set)
    P.Patch.TD(3:20) = repmat(pars(123),1,18);
end
% Ls0Pas (19)
if contains('Ls0Pas_Lv',par_set)
    P.Patch.Ls0Pas(3:20) = pars(141:158);
end
if contains('Ls0Pas_g_Lv',par_set)
    P.Patch.Ls0Pas(3:20) = repmat(pars(141),1,18);
end
P.Patch.Ls0Pas(1) = pars(159);
% SfPas (18)
if contains('SfPas_Lv',par_set)
    P.Patch.SfPas(3:20) = pars(160:177);
end
if contains('SfPas_g_Lv',par_set)
    P.Patch.SfPas(3:14) = repmat(pars(160)*2.2311e4,1,12);
    P.Patch.SfPas(15:20) = repmat(pars(160)*2.2086e4,1,6);
end
% k1 (19)
if contains('k1_Lv',par_set)
    P.Patch.k1(3:20) = pars(178:195);
end
if contains('k1_g_Lv',par_set)
    P.Patch.k1(3:20) = repmat(pars(178),1,18);
end
P.Patch.k1(1) = pars(196);
% dLsPas (18)
if contains('dLsPas_Lv',par_set)
    P.Patch.dLsPas(3:20) = pars(197:214);
end
if contains('dLsPas_g_Lv',par_set)
    P.Patch.dLsPas(3:20) = repmat(pars(197),1,18);
end
% LenSeriesElement (1)
if contains('LenSeriesElement_g_Lv',par_set)
    P.Patch.LenSeriesElement(3:20) = repmat(pars(215),1,18);
end
% ADO (19)
if contains('ADO_Lv',par_set)
    P.Patch.ADO(3:20) = pars(216:233);
end
if contains('ADO_g_Lv',par_set)
    P.Patch.ADO(3:20) = repmat(pars(216),1,18);
end
P.Patch.ADO(1) = pars(234);
% LDAD (18)
if contains('LDAD_Lv',par_set)
    P.Patch.LDAD(3:20) = pars(235:252);
end
if contains('LDAD_g_Lv',par_set)
    P.Patch.LDAD(3:20) = repmat(pars(235),1,18);
end
% LDCI (18)
if contains('LDCI_Lv',par_set)
    P.Patch.LDCI(3:20) = pars(253:270);
end
if contains('LDCI_g_Lv',par_set)
    P.Patch.LDCI(3:20) = repmat(pars(253),1,18);
end