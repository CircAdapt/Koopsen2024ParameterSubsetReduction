% ---MorrisScreening_input.m---
% 
% Defines all input of the Morris Screening Method, after which
% trajectories are evaluated by MorrisScreening_analysis.m
%
% Koopsen T.
% Last modified: 11/27/2023

% Add path with CircAdapt C++ code
addpath([pwd,'\CircAdapt-master\matlab\example']);

% Load reference P-struct
load('PRefHF.mat')

% Define reference wall volumes of left ventricle (LV), 
% right ventricle (RV), left atrium (LA), and right atrium (RA)
VWall_LV_ref = sum(P.Patch.VWall(3:4));
VWall_RV_ref = P.Patch.VWall(5);
VWall_LA_ref = P.Patch.VWall(1);
VWall_RA_ref = P.Patch.VWall(2);
% Define reference wall areas of left ventricle (LV), 
% right ventricle (RV), left atrium (LA), right atrium (RA)
AmRef_LV_ref = sum(P.Patch.AmRef(3:4));
AmRef_RV_ref = P.Patch.AmRef(5);
AmRef_LA_ref = P.Patch.AmRef(1);
AmRef_RA_ref = P.Patch.AmRef(2);

% Subdivide LV wall into 6 patches (6-segment model)
% Names used in this script:                              Sv1 - AntSept, Sv2 - Sept, Lv1 - Ant, Lv2 - Inf, Lv3 - Lat, Lv4 - Post
% Names used in Koopsen et al. (2024) Biomed. Eng. Onl.:  S1  - AntSept, S2  - Sept, LV1 - Ant, LV2 - Inf, LV3 - Lat, LV4 - Post
SplitMerge('Sv1',2);
SplitMerge('Lv1',4);

% Run the 6-segment model
RunCPP

% Define all parameter limits (npar=174) for the Morris Screening Method
all_lims = zeros(174,4);

% Global circulatory parameters (P.General) (1-3)
q0_lim = [2.1,5.5]/(60*1000); % limits of study population (n=13)
all_lims(1,1:2) = q0_lim;
p0_lim = P.General.p0*[0.75,1.25];
all_lims(2,1:2) = p0_lim;
tCycle_lim = 60./[91.1,48.4]; % limits of study population (n=13)
all_lims(3,1:2) = tCycle_lim;

% ArtVen parameters (4-15)
k_Sy_lim = GetFt('ArtVen','k','Sy').*[0.8,1.2];
all_lims(4,:) = [k_Sy_lim(1),k_Sy_lim(3),k_Sy_lim(2),k_Sy_lim(4)];
k_Pu_lim = GetFt('ArtVen','k','Pu').*[0.8,1.2];
all_lims(5,:) = [k_Pu_lim(1),k_Pu_lim(3),k_Pu_lim(2),k_Pu_lim(4)];
Len_Sy_lim = GetFt('ArtVen','Len','Sy').*[0.8,1.2];
all_lims(6,:) = [Len_Sy_lim(1),Len_Sy_lim(3),Len_Sy_lim(2),Len_Sy_lim(4)];
Len_Pu_lim = GetFt('ArtVen','Len','Pu').*[0.8,1.2];
all_lims(7,:) = [Len_Pu_lim(1),Len_Pu_lim(3),Len_Pu_lim(2),Len_Pu_lim(4)];
p0_Sy_lim = GetFt('ArtVen','p0','Sy').*[0.8,1.2];
all_lims(8,:) = [p0_Sy_lim(1),p0_Sy_lim(3),p0_Sy_lim(2),p0_Sy_lim(4)];
p0_Pu_lim = GetFt('ArtVen','p0','Pu').*[0.8,1.2];
all_lims(9,:) = [p0_Pu_lim(1),p0_Pu_lim(3),p0_Pu_lim(2),p0_Pu_lim(4)];
A0_Sy_lim = GetFt('ArtVen','A0','Sy').*[0.8,1.2];
all_lims(10,:) = [A0_Sy_lim(1),A0_Sy_lim(3),A0_Sy_lim(2),A0_Sy_lim(4)];
A0_Pu_lim = GetFt('ArtVen','A0','Pu').*[0.8,1.2];
all_lims(11,:) = [A0_Pu_lim(1),A0_Pu_lim(3),A0_Pu_lim(2),A0_Pu_lim(4)];
AWall_Sy_lim = GetFt('ArtVen','AWall','Sy').*[0.8,1.2];
all_lims(12,:) = [AWall_Sy_lim(1),AWall_Sy_lim(3),AWall_Sy_lim(2),AWall_Sy_lim(4)];
AWall_Pu_lim = GetFt('ArtVen','AWall','Pu').*[0.8,1.2];
all_lims(13,:) = [AWall_Pu_lim(1),AWall_Pu_lim(3),AWall_Pu_lim(2),AWall_Pu_lim(4)];
kAV_Sy_lim = GetFt('ArtVen','kAV','Sy')*[0.8,1.2];
all_lims(14,1:2) = kAV_Sy_lim;
kAV_Pu_lim = GetFt('ArtVen','kAV','Pu')*[0.8,1.2];
all_lims(15,1:2) = kAV_Pu_lim;

% Pericardial parameters (16-18)
k_Peri_lim = GetFt('Bag','k','Peri')*[0.8,1.2];
all_lims(16,1:2) = k_Peri_lim;
VRef_Peri_lim = GetFt('Bag','VRef','Peri')*[0.8,1.2];
all_lims(17,1:2) = VRef_Peri_lim;
pAdapt_Peri_lim = GetFt('Bag','pAdapt','Peri')*[0.8,1.2];
all_lims(18,1:2) = pAdapt_Peri_lim;

% Valvular parameters (19-26)
AOpen_RaRv_lim = GetFt('Valve','AOpen','RaRv')*[0.8,1.2];
all_lims(19,1:2) = AOpen_RaRv_lim;
AOpen_RvPuArt_lim = GetFt('Valve','AOpen','RvPuArt')*[0.8,1.2];
all_lims(20,1:2) = AOpen_RvPuArt_lim;
AOpen_LaLv_lim = GetFt('Valve','AOpen','LaLv')*[0.8,1.2];
all_lims(21,1:2) = AOpen_LaLv_lim;
AOpen_LvSyArt_lim = GetFt('Valve','AOpen','LvSyArt')*[0.8,1.2];
all_lims(22,1:2) = AOpen_LvSyArt_lim;

Len_RaRv_lim = GetFt('Valve','Len','RaRv')*[0.8,1.2];
all_lims(23,1:2) = Len_RaRv_lim;
Len_RvPuArt_lim = GetFt('Valve','Len','RvPuArt')*[0.8,1.2];
all_lims(24,1:2) = Len_RvPuArt_lim;
Len_LaLv_lim = GetFt('Valve','Len','LaLv')*[0.8,1.2];
all_lims(25,1:2) = Len_LaLv_lim;
Len_LvSyArt_lim = GetFt('Valve','Len','LvSyArt')*[0.8,1.2];
all_lims(26,1:2) = Len_LvSyArt_lim;

% Electrical substrate (dT) - 8 parameters (27-34)
% 
% For parameter definitions see Koopsen et al. (2024) Biomed.Eng.Onl. Figure S3

dTauAv_lim = [-0.050,0.100];
all_lims(27,1:2) = dTauAv_lim;
tau1_lim = [-0.030,0.030]; % interventricular delay (RV->septum)
all_lims(28,1:2) = tau1_lim;
tau2_lim = [0,0.120]; % intraventricular delay (septum->LV free wall)
all_lims(29,1:2) = tau2_lim;
alfa1_lim = [1/6,3/6];
all_lims(30,1:2) = alfa1_lim;
alfa2_lim = [3/6,5/6];
all_lims(31,1:2) = alfa2_lim;
beta1_lim = [1/6,3/6];
all_lims(32,1:2) = beta1_lim;
beta2_lim = [3/6,5/6];
all_lims(33,1:2) = beta2_lim;

dT_La_lim = [-0.010,0.050]; % left atrium
all_lims(34,1:2) = dT_La_lim;

% Patch parameters - 14 parameters*10 (10 = global offset (GO) of LV, 6 LV patches (p),RV,RA,LA)

% VWall
VWall_LV_GO_lim = [0.8,1.2]*VWall_LV_ref;
all_lims(35,1:2) = VWall_LV_GO_lim;
VWall_LV_p_mf_lim = [0.8,1.2];
all_lims(36:41,1:2) = repmat(VWall_LV_p_mf_lim,6,1);
VWall_RV_lim = [0.8,1.2]*VWall_RV_ref;
all_lims(42,1:2) = VWall_RV_lim;
VWall_LA_lim = [0.8,1.2]*VWall_LA_ref;
all_lims(43,1:2) = VWall_LA_lim;
VWall_RA_lim = [0.8,1.2]*VWall_RA_ref;
all_lims(44,1:2) = VWall_RA_lim;

% AmRef
AmRef_LV_GO_lim = [0.8,1.2]*AmRef_LV_ref;
all_lims(45,1:2) = AmRef_LV_GO_lim;
AmRef_LV_p_mf_lim = [0.8,1.2];
all_lims(46:51,1:2) = repmat(AmRef_LV_p_mf_lim,6,1);
AmRef_RV_lim = [0.8,1.2]*AmRef_RV_ref;
all_lims(52,1:2) = AmRef_RV_lim;
AmRef_LA_lim = [0.8,1.2]*AmRef_LA_ref;
all_lims(53,1:2) = AmRef_LA_lim;
AmRef_RA_lim = [0.8,1.2]*AmRef_RA_ref;
all_lims(54,1:2) = AmRef_RA_lim;

% SfAct
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
SfAct_LV_GO_lim = [0.5,1.2]*GetFt('Patch','SfAct','Lv1');
all_lims(55,1:2) = SfAct_LV_GO_lim;
SfAct_LV_p_mf_I_lim = [0,1.2];
all_lims([56,58,60],1:2) = repmat(SfAct_LV_p_mf_I_lim,3,1);
SfAct_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([57,59,61],1:2) = repmat(SfAct_LV_p_mf_NI_lim,3,1);
SfAct_RV_lim = GetFt('Patch','SfAct','Rv1')*[0.5,1.2];
all_lims(62,1:2) = SfAct_RV_lim;
SfAct_a_lim = GetFt('Patch','SfAct','La1')*[0.5,1.2];
all_lims(63:64,1:2) = repmat(SfAct_a_lim,2,1);

% vMax
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
vMax_LV_GO_lim = [0.8,1.2]*GetFt('Patch','vMax','Lv1');
all_lims(65,1:2) = vMax_LV_GO_lim;
vMax_LV_p_mf_I_lim = [0.5,1.2];
all_lims([66,68,70],1:2) = repmat(vMax_LV_p_mf_I_lim,3,1);
vMax_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([67,69,71],1:2) = repmat(vMax_LV_p_mf_NI_lim,3,1);
vMax_RV_lim = GetFt('Patch','vMax','Rv1')*[0.8,1.2];
all_lims(72,1:2) = vMax_RV_lim;
vMax_a_lim = GetFt('Patch','vMax','La1')*[0.8,1.2];
all_lims(73:74,1:2) = repmat(vMax_a_lim,2,1);

% TR
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
TR_LV_GO_lim = [0.8,1.2]*GetFt('Patch','TR','Lv1');
all_lims(75,1:2) = TR_LV_GO_lim;
TR_LV_p_mf_I_lim = [0.8,1.2];
all_lims([76,78,80],1:2) = repmat(TR_LV_p_mf_I_lim,3,1);
TR_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([77,79,81],1:2) = repmat(TR_LV_p_mf_NI_lim,3,1);
TR_RV_lim = GetFt('Patch','TR','Rv1')*[0.8,1.2];
all_lims(82,1:2) = TR_RV_lim;
TR_a_lim = GetFt('Patch','TR','La1')*[0.8,1.2];
all_lims(83:84,1:2) = repmat(TR_a_lim,2,1);

% TD
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
TD_LV_GO_lim = [0.8,1.2]*GetFt('Patch','TD','Lv1');
all_lims(85,1:2) = TD_LV_GO_lim;
TD_LV_p_mf_I_lim = [0.8,1.2];
all_lims([86,88,90],1:2) = repmat(TD_LV_p_mf_I_lim,3,1);
TD_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([87,89,91],1:2) = repmat(TD_LV_p_mf_NI_lim,3,1);
TD_RV_lim = GetFt('Patch','TD','Rv1')*[0.8,1.2];
all_lims(92,1:2) = TD_RV_lim;
TD_a_lim = GetFt('Patch','TD','La1')*[0.8,1.2];
all_lims(93:94,1:2) = repmat(TD_a_lim,2,1);

% Ls0Pas
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
Ls0Pas_LV_GO_lim = [0.95,1.05]*GetFt('Patch','Ls0Pas','Lv1');
all_lims(95,1:2) = Ls0Pas_LV_GO_lim;
Ls0Pas_LV_p_mf_I_lim = [0.95,1.05];
all_lims([96,98,100],1:2) = repmat(Ls0Pas_LV_p_mf_I_lim,3,1);
Ls0Pas_LV_p_mf_NI_lim = [0.95,1.05];
all_lims([97,99,101],1:2) = repmat(Ls0Pas_LV_p_mf_NI_lim,3,1);
Ls0Pas_RV_lim = GetFt('Patch','Ls0Pas','Rv1')*[0.95,1.05];
all_lims(102,1:2) = Ls0Pas_RV_lim;
Ls0Pas_a_lim = GetFt('Patch','Ls0Pas','La1')*[0.95,1.05];
all_lims(103:104,1:2) = repmat(Ls0Pas_a_lim,2,1);

% SfPas
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
SfPas_LV_GO_mf_lim = [0.8,1.2];
all_lims(105,1:2) = SfPas_LV_GO_mf_lim;
SfPas_LV_p_mf_I_lim = [0.8,10];
all_lims([106,108,110],1:2) = repmat(SfPas_LV_p_mf_I_lim,3,1);
SfPas_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([107,109,111],1:2) = repmat(SfPas_LV_p_mf_NI_lim,3,1);
SfPas_RV_lim = GetFt('Patch','SfPas','Rv1')*[0.8,1.2];
all_lims(112,1:2) = SfPas_RV_lim;
SfPas_a_lim = GetFt('Patch','SfPas','La1')*[0.8,1.2];
all_lims(113:114,1:2) = repmat(SfPas_a_lim,2,1);

% k1
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
k1_LV_GO_lim = [0.8,1.2]*GetFt('Patch','k1','Lv1');
all_lims(115,1:2) = k1_LV_GO_lim;
k1_LV_p_mf_I_lim = [0.8,3];
all_lims([116,118,120],1:2) = repmat(k1_LV_p_mf_I_lim,3,1);
k1_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([117,119,121],1:2) = repmat(k1_LV_p_mf_NI_lim,3,1);
k1_RV_lim = [8,12];
all_lims(122,1:2) = k1_RV_lim;
k1_a_lim = [8,12];
all_lims(123:124,1:2) = repmat(k1_a_lim,2,1);

% dLsPas
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
dLsPas_LV_GO_lim = [0.8,1.2]*GetFt('Patch','dLsPas','Lv1');
all_lims(125,1:2) = dLsPas_LV_GO_lim;
dLsPas_LV_p_mf_I_lim = [0.8,1.2];
all_lims([126,128,130],1:2) = repmat(dLsPas_LV_p_mf_I_lim,3,1);
dLsPas_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([127,129,131],1:2) = repmat(dLsPas_LV_p_mf_NI_lim,3,1);
dLsPas_RV_lim = [0.8,1.2]*GetFt('Patch','dLsPas','Rv1');
all_lims(132,1:2) = dLsPas_RV_lim;
dLsPas_a_lim = [0.8,1.2]*GetFt('Patch','dLsPas','La1');
all_lims(133:134,1:2) = repmat(dLsPas_a_lim,2,1);

% LenSeriesElement (LenSE)
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
LenSE_LV_GO_lim = [0.8,1.2]*GetFt('Patch','LenSeriesElement','Lv1');
all_lims(135,1:2) = LenSE_LV_GO_lim;
LenSE_LV_p_mf_I_lim = [0.8,1.2];
all_lims([136,138,140],1:2) = repmat(LenSE_LV_p_mf_I_lim,3,1);
LenSE_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([137,139,141],1:2) = repmat(LenSE_LV_p_mf_NI_lim,3,1);
LenSE_RV_lim = [0.8,1.2]*GetFt('Patch','LenSeriesElement','Rv1');
all_lims(142,1:2) = LenSE_RV_lim;
LenSE_a_lim = [0.8,1.2]*GetFt('Patch','LenSeriesElement','La1');
all_lims(143:144,1:2) = repmat(LenSE_a_lim,2,1);

% Activation duration offset (ADO)
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
ADO_LV_GO_lim = [0.8,1.2]*GetFt('Patch','ADO','Lv1');
all_lims(145,1:2) = ADO_LV_GO_lim;
ADO_LV_p_mf_I_lim = [0.8,1.2];
all_lims([146,148,150],1:2) = repmat(ADO_LV_p_mf_I_lim,3,1);
ADO_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([147,149,151],1:2) = repmat(ADO_LV_p_mf_NI_lim,3,1);
ADO_RV_lim = [0.8,1.2]*GetFt('Patch','ADO','Rv1');
all_lims(152,1:2) = ADO_RV_lim;
ADO_a_lim = [0.8,1.2]*GetFt('Patch','ADO','La1');
all_lims(153:154,1:2) = repmat(ADO_a_lim,2,1);

% Length dependency of activation duration (LDAD)
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
LDAD_LV_GO_lim = [0.8,1.2]*GetFt('Patch','LDAD','Lv1');
all_lims(155,1:2) = LDAD_LV_GO_lim;
LDAD_LV_p_mf_I_lim = [0.5,1.2];
all_lims([156,158,160],1:2) = repmat(LDAD_LV_p_mf_I_lim,3,1);
LDAD_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([157,159,161],1:2) = repmat(LDAD_LV_p_mf_NI_lim,3,1);
LDAD_RV_lim = [0.8,1.2]*GetFt('Patch','LDAD','Rv1');
all_lims(162,1:2) = LDAD_RV_lim;
LDAD_a_lim = [0.8,1.2]*GetFt('Patch','LDAD','La1');
all_lims(163:164,1:2) = repmat(LDAD_a_lim,2,1);

% Length dependency of contractility increase (LDCI)
% I = 'ischemic' patches (S1,LV1,LV3)
% NI = 'non-ischemic' patches (S2,LV2,LV4)
LDCI_LV_GO_lim = [0.8,1.2]*GetFt('Patch','LDCI','Lv1');
all_lims(165,1:2) = LDCI_LV_GO_lim;
LDCI_LV_p_mf_I_lim = [0.5,1.2];
all_lims([166,168,170],1:2) = repmat(LDCI_LV_p_mf_I_lim,3,1);
LDCI_LV_p_mf_NI_lim = [0.8,1.2];
all_lims([167,169,171],1:2) = repmat(LDCI_LV_p_mf_NI_lim,3,1);
LDCI_RV_lim = [0.8,1.2]*GetFt('Patch','LDCI','Rv1');
all_lims(172,1:2) = LDCI_RV_lim;
LDCI_a_lim = [0.8,1.2]*GetFt('Patch','LDCI','La1');
all_lims(173:174,1:2) = repmat(LDCI_a_lim,2,1);

fields1 = [repmat({'General'},3,1);repmat({'ArtVen'},12,1);repmat({'Bag'},3,1);repmat({'Valve'},8,1);'General';repmat({'Patch'},147,1)];
fields2 = ['q0';'p0';'tCycle';repmat({'k'},2,1);repmat({'Len'},2,1);repmat({'p0'},2,1);repmat({'A0'},2,1);repmat({'AWall'},2,1);repmat({'kAV'},2,1);'k';'VRef';'pAdapt';repmat({'AOpen'},4,1);repmat({'Len'},4,1);'dTauAv';'tau1';'tau2';repmat({'dTtau1'},4,1);'dT';repmat({'VWall'},10,1);repmat({'AmRef'},10,1);repmat({'SfAct'},10,1);repmat({'vMax'},10,1);repmat({'TR'},10,1);repmat({'TD'},10,1);repmat({'Ls0Pas'},10,1);repmat({'SfPas'},10,1);repmat({'k1'},10,1);repmat({'dLsPas'},10,1);repmat({'LenSeriesElement'},10,1);repmat({'ADO'},10,1);repmat({'LDAD'},10,1);repmat({'LDCI'},10,1)];
fields3 = [repmat({'glob'},3,1);repmat({'Sy';'Pu'},6,1);repmat({'Peri'},3,1);repmat({'RaRv';'RvPuArt';'LaLv';'LvSyArt'},2,1);'glob';'E';'L';'ME';'ML';'ME';'ML';'La1';repmat({'GO';'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'Rv1';'La1';'Ra1'},14,1)];

% Simulate trajectories
nsubsteps_sp = 2; % number of substeps from reference to starting point (to enhance stability, but increases computational time)
nsubsteps_tr = 5; % number of substeps within the trajectory itself (to enhance stability, but increases computational time)

ntrajectories = 30000; % total number of trajectories which are evaluated

% Define the random trajectories
starting_points_0 = (randi(8,ntrajectories,174)-1)/7;
starting_points_0 = starting_points_0';
trajectories = zeros(174,ntrajectories);
for i = 1:ntrajectories
    trajectories(:,i) = randperm(174,174)';
end

% Enter output folder
savefolder = [pwd,'\Morris screening output'];

% Run all trajectories in parallel (function: MorrisScreening_analysis.m)
parfor n = 1:ntrajectories
    [delta_strain_indices_LV,delta_LVEDV,delta_LVSV] = MorrisScreening_analysis('PRefHF.mat',all_lims,fields1,fields2,fields3,nsubsteps_sp,nsubsteps_tr,starting_points_0(:,n),trajectories(:,n),savefolder,n);
end


