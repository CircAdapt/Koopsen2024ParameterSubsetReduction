% ---setBounds.m---
%
% Set parameter boundaries for Monte Carlo simulations
%
% Koopsen T.
% Last modified: 11/27/2023

function par_bounds = setBounds(refP,par_set,measurement)

% Input arguments:
% refP        - reference P-struct
% par_set     - parameter set used for optimization
% measurement - struct containing patient-specific input data
%
% Output arguments:
% par_bounds  - matrix(270*2) with parameter boundaries, row indexing as follows:
% 1            - q0
% 2            - p0
% 3            - k_Sy
% 4            - p0_Sy
% 5            - A0_Sy
% 6            - dTauAv
% 7-24         - dT_Lv (18 patches)
% 25           - dT_La
% 26-43        - VWall_Lv (18 patches)
% 44           - VWall_La
% 45-62        - AmRef_Lv (18 patches)
% 63           - AmRef_Rv
% 64           - AmRef_La
% 65-82        - SfAct_Lv (18 patches)
% 83           - SfAct_Rv
% 84           - SfAct_La
% 85-102       - vMax_Lv (18 patches)
% 103-120      - TR_Lv (18 patches)
% 121          - TR_Rv
% 122          - TR_La
% 123-140      - TD_Lv (18 patches)
% 141-158      - Ls0Pas_Lv (18 patches)
% 159          - Ls0Pas_La
% 160-177      - SfPas_Lv (18 patches)
% 178-195      - k1_Lv (18 patches)
% 196          - k1_La
% 197-214      - dLsPas_Lv (18 patches)
% 215          - LenSE_Lv_G (global)
% 216-233      - ADO_Lv (18 patches)
% 234          - ADO_La
% 235-252      - LDAD_Lv (18 patches)
% 253-270      - LDCI_Lv (18 patches)

% Load reference P-struct
load(refP,'P')

par_bounds = zeros(270,2); % 1st column = lower boundary, 2nd column = upper boundary

% Global circulation (2)
% q0 from measurement (q0_mea)
cycle_duration = measurement.time(end);
SV_mea = measurement.LV.EDV-measurement.LV.ESV;
q0_mea = ((60/cycle_duration)*SV_mea)/(60*1000*1000);
if contains('q0',par_set)
    par_bounds(1,1:2) = q0_mea*[0.8,1.2];
else
    par_bounds(1,1:2) = q0_mea*[1,1];
end
% p0
if contains('p0',par_set)
    par_bounds(2,1:2) = P.General.p0*[0.75,1.25];
else
    par_bounds(2,1:2) = P.General.p0*[1,1];
end

% ArtVen (3)
% k_Sy
if contains('k_Sy',par_set)
    par_bounds(3,:) = [0.8,1.2];
else
    par_bounds(3,:) = [1,1];
end
% p0_Sy
if contains('p0_Sy',par_set)
    par_bounds(4,:) = [0.8,1.2];
else
    par_bounds(4,:) = [1,1];
end
% A0_Sy
if contains('A0_Sy',par_set)
    par_bounds(5,:) = [0.8,1.2];
else
    par_bounds(5,:) = [1,1];
end

% Electrical substrate (20)
% dTauAv
if contains('dTauAv',par_set)
    par_bounds(6,1:2) = [-0.050,0.130];
else
    par_bounds(6,1:2) = [0,0];
end

% dT_Lv (18 patches)
if contains('dT_Lv',par_set)
    par_bounds(7:18,1:2) = repmat([0,0.120],12,1); % LVfw
    par_bounds(19:24,1:2) = repmat([-0.030,0.060],6,1); % septum
elseif contains('dT_g_Lv',par_set) % global dT
    par_bounds(7,1:2) = [0,0.120];
else
    par_bounds(7:18,1:2) = repmat([0,0],12,1); % LVfw
    par_bounds(19:24,1:2) = repmat([0,0],6,1); % septum
end

% dT_La
if contains('dT_La',par_set)
    par_bounds(25,1:2) = [-0.010,0.050];
else
    par_bounds(25,1:2) = [0,0];
end

% Other patch parameters (245)
% VWall_Lv (18 patches)
if contains('VWall_Lv',par_set)
    par_bounds(26:37,1:2) = repmat(P.Patch.VWall(3)*[0.64,1.44],12,1); % LVfw
    par_bounds(38:43,1:2) = repmat(P.Patch.VWall(15)*[0.64,1.44],6,1); % septum
% VWall_Lv_g (global)
elseif contains('VWall_g_Lv',par_set)
    par_bounds(26,1:2) = [0.8,1.2];
    par_bounds(27:43,1:2) = repmat([1,1],17,1);
else
    par_bounds(26:37,1:2) = repmat(P.Patch.VWall(3)*[1,1],12,1); % LVfw
    par_bounds(38:43,1:2) = repmat(P.Patch.VWall(15)*[1,1],6,1); % septum
end

% VWall_La
if contains('VWall_La',par_set)
    par_bounds(44,1:2) = P.Patch.VWall(1)*[0.8,1.2];
else
    par_bounds(44,1:2) = P.Patch.VWall(1)*[1,1];
end

% AmRef_Lv (18 patches)
if contains('AmRef_Lv',par_set)
    par_bounds(45:56,1:2) = repmat(P.Patch.AmRef(3)*[0.64,1.44],12,1); % LVfw
    par_bounds(57:62,1:2) = repmat(P.Patch.AmRef(15)*[0.64,1.44],6,1); % septum
% AmRef_Lv_g (global)
elseif contains('AmRef_g_Lv',par_set)
    par_bounds(45,1:2) = [0.8,1.2];
    par_bounds(46:62,1:2) = repmat([1,1],17,1);
else
    par_bounds(45:56,1:2) = repmat(P.Patch.AmRef(3)*[1,1],12,1); % LVfw
    par_bounds(57:62,1:2) = repmat(P.Patch.AmRef(15)*[1,1],6,1); % septum
end

% AmRef_Rv
if contains('AmRef_Rv',par_set)
    par_bounds(63,1:2) = P.Patch.AmRef(21)*[0.8,1.2];
else
    par_bounds(63,1:2) = P.Patch.AmRef(21)*[1,1];
end

% AmRef_La
if contains('AmRef_La',par_set)
    par_bounds(64,1:2) = P.Patch.AmRef(1)*[0.8,1.2];
else
    par_bounds(64,1:2) = P.Patch.AmRef(1)*[1,1];
end

% SfAct_Lv (18 patches)
if contains('SfAct_Lv',par_set)
    par_bounds(65:82,1:2) = repmat(P.Patch.SfAct(3)*[0,1.44],18,1);
% SfAct_Lv_g (global)
elseif contains('SfAct_g_Lv',par_set)
    par_bounds(65:82,1:2) = repmat(P.Patch.SfAct(3)*[0.5,1.2],18,1);
    par_bounds(66:82,1:2) = repmat(P.Patch.SfAct(3)*[1,1],17,1);
else
    par_bounds(65:82,1:2) = repmat(P.Patch.SfAct(3)*[1,1],18,1);
end

% SfAct_Rv
if contains('SfAct_Rv',par_set)
    par_bounds(83,1:2) = P.Patch.SfAct(21)*[0.5,1.2];
else
    par_bounds(83,1:2) = P.Patch.SfAct(21)*[1,1];
end

% SfAct_La
if contains('SfAct_La',par_set)
    par_bounds(84,1:2) = P.Patch.SfAct(1)*[0.5,1.2];
else
    par_bounds(84,1:2) = P.Patch.SfAct(1)*[1,1];
end

% vMax_Lv (18 patches)
if contains('vMax_Lv',par_set)
    par_bounds(85:102,1:2) = repmat(P.Patch.vMax(3)*[0.4,1.44],18,1);
elseif contains('vMax_g_Lv',par_set)
    par_bounds(85:102,1:2) = repmat(P.Patch.vMax(3)*[0.8,1.2],18,1);
    par_bounds(86:102,1:2) = repmat(P.Patch.vMax(3)*[1,1],17,1);
else
    par_bounds(85:102,1:2) = repmat(P.Patch.vMax(3)*[1,1],18,1);
end

% TR_Lv (18 patches)
if contains('TR_Lv',par_set)
    par_bounds(103:120,1:2) = repmat(P.Patch.TR(3)*[0.64,1.44],18,1);
elseif contains('TR_g_Lv',par_set)
    par_bounds(103:120,1:2) = repmat(P.Patch.TR(3)*[0.8,1.2],18,1);
    par_bounds(104:120,1:2) = repmat(P.Patch.TR(3)*[1,1],17,1);
else
    par_bounds(103:120,1:2) = repmat(P.Patch.TR(3)*[1,1],18,1);
end

% TR_Rv
if contains('TR_Rv',par_set)
    par_bounds(121,1:2) = P.Patch.TR(21)*[0.8,1.2];
else
    par_bounds(121,1:2) = P.Patch.TR(21)*[1,1];
end

% TR_La
if contains('TR_La',par_set)
    par_bounds(122,1:2) = P.Patch.TR(1)*[0.8,1.2];
else
    par_bounds(122,1:2) = P.Patch.TR(1)*[1,1];
end

% TD_Lv (18 patches)
if contains('TD_Lv',par_set)
    par_bounds(123:140,1:2) = repmat(P.Patch.TD(3)*[0.64,1.44],18,1);
elseif contains('TD_g_Lv',par_set)
    par_bounds(123:140,1:2) = repmat(P.Patch.TD(3)*[0.8,1.2],18,1);
    par_bounds(124:140,1:2) = repmat(P.Patch.TD(3)*[1,1],17,1);
else
    par_bounds(123:140,1:2) = repmat(P.Patch.TD(3)*[1,1],18,1);
end

% Ls0Pas_Lv (18 patches)
if contains('Ls0Pas_Lv',par_set)
    par_bounds(141:158,1:2) = repmat(P.Patch.Ls0Pas(3)*[0.9025,1.1025],18,1);
elseif contains('Ls0Pas_g_Lv',par_set)
    par_bounds(141:158,1:2) = repmat(P.Patch.Ls0Pas(3)*[0.95,1.05],18,1);
    par_bounds(142:158,1:2) = repmat(P.Patch.Ls0Pas(3)*[1,1],17,1);
else
    par_bounds(141:158,1:2) = repmat(P.Patch.Ls0Pas(3)*[1,1],18,1);
end

% Ls0Pas_La
if contains('Ls0Pas_La',par_set)
    par_bounds(159,1:2) = P.Patch.Ls0Pas(1)*[0.95,1.05];
else
    par_bounds(159,1:2) = P.Patch.Ls0Pas(1)*[1,1];
end

% SfPas_Lv (18 patches)
if contains('SfPas_Lv',par_set)
    par_bounds(160:171,1:2) = repmat(P.Patch.SfPas(3)*[0.64,12],12,1); % LVfw
    par_bounds(172:177,1:2) = repmat(P.Patch.SfPas(15)*[0.64,12],6,1); % septum
elseif contains('SfPas_g_Lv',par_set)
    par_bounds(160,1:2) = [0.8,1.2];
    par_bounds(161:177,1:2) = repmat([1,1],17,1);
else
    par_bounds(160:171,1:2) = repmat(P.Patch.SfPas(3)*[1,1],12,1); % LVfw
    par_bounds(172:177,1:2) = repmat(P.Patch.SfPas(15)*[1,1],6,1); % septum
end

% k1_Lv (18 patches)
if contains('k1_Lv',par_set)
    par_bounds(178:195,1:2) = repmat(P.Patch.k1(3)*[0.64,3.6],18,1);
elseif contains('k1_g_Lv',par_set)
    par_bounds(178:195,1:2) = repmat(P.Patch.k1(3)*[0.8,1.2],18,1);
    par_bounds(179:195,1:2) = repmat(P.Patch.k1(3)*[1,1],17,1);
else
    par_bounds(178:195,1:2) = repmat(P.Patch.k1(3)*[1,1],18,1);
end

% k1_La
if contains('k1_La',par_set)
    par_bounds(196,1:2) = P.Patch.k1(1)*[0.8,1.2];
else
    par_bounds(196,1:2) = P.Patch.k1(1)*[1,1];
end

% dLsPas_Lv (18 patches)
if contains('dLsPas_Lv',par_set)
    par_bounds(197:214,1:2) = repmat(P.Patch.dLsPas(3)*[0.64,1.44],18,1);
elseif contains('dLsPas_g_Lv',par_set)
    par_bounds(197:214,1:2) = repmat(P.Patch.dLsPas(3)*[0.8,1.2],18,1);
    par_bounds(198:214,1:2) = repmat(P.Patch.dLsPas(3)*[1,1],17,1);
else
    par_bounds(197:214,1:2) = repmat(P.Patch.dLsPas(3)*[1,1],18,1);
end

% LenSE_Lv_g (global)
if contains('LenSeriesElement_g_Lv',par_set)
    par_bounds(215,1:2) = P.Patch.LenSeriesElement(3)*[0.8,1.2];
else
    par_bounds(215,1:2) = P.Patch.LenSeriesElement(3)*[1,1];
end

% ADO_Lv (18 patches)
if contains('ADO_Lv',par_set)
    par_bounds(216:233,1:2) = repmat(P.Patch.ADO(3)*[0.64,1.44],18,1);
elseif contains('ADO_g_Lv',par_set)
    par_bounds(216:233,1:2) = repmat(P.Patch.ADO(3)*[0.8,1.2],18,1);
    par_bounds(217:233,1:2) = repmat(P.Patch.ADO(3)*[1,1],17,1);
else
    par_bounds(216:233,1:2) = repmat(P.Patch.ADO(3)*[1,1],18,1);
end

% ADO_La
if contains('ADO_La',par_set)
    par_bounds(234,1:2) = P.Patch.ADO(1)*[0.8,1.2];
else
    par_bounds(234,1:2) = P.Patch.ADO(1)*[1,1];
end

% LDAD_Lv (18 patches)
if contains('LDAD_Lv',par_set)
    par_bounds(235:252,1:2) = repmat(P.Patch.LDAD(3)*[0.4,1.44],18,1);
elseif contains('LDAD_g_Lv',par_set)
    par_bounds(235:252,1:2) = repmat(P.Patch.LDAD(3)*[0.8,1.2],18,1);
    par_bounds(236:252,1:2) = repmat(P.Patch.LDAD(3)*[1,1],17,1);
else
    par_bounds(235:252,1:2) = repmat(P.Patch.LDAD(3)*[1,1],18,1);
end

% LDCI_Lv (18 patches)
if contains('LDCI_Lv',par_set)
    par_bounds(253:270,1:2) = repmat(P.Patch.LDCI(3)*[0.4,1.44],18,1);
elseif contains('LDCI_g_Lv',par_set)
    par_bounds(253:270,1:2) = repmat(P.Patch.LDCI(3)*[0.8,1.2],18,1);
    par_bounds(254:270,1:2) = repmat(P.Patch.LDCI(3)*[1,1],17,1);
else
    par_bounds(253:270,1:2) = repmat(P.Patch.LDCI(3)*[1,1],18,1);
end

% tCycle
% if contains('tCycle',par_set)
%     %par_bounds(271,1:2) = 60./[91.1,48.4];
%     par_bounds(271,1:2) = measurement.cdmean + measurement.cdstd*[-2,2];
% else
%     par_bounds(271,1:2) = [0.85,0.85];
% end