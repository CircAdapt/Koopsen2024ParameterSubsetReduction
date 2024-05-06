% ---SobolSamples.m---
%
% Initialize Sobol samples, which are then evaluated by RunSobolSimulation.m
%
% Koopsen T.
% Last modified: 11/27/2023

% Parameters (103) are based on the 6-segment LV model:
% - Global hemodynamics (3)
%   1 - q0
%   2 - p0
%   3 - tCycle
% - ArtVen (3)
%   4 - kSy
%   5 - p0Sy
%   6 - A0Sy
% - Patch (97) 
%   - Electrical substrate (8)
%     7 - dTauAv
%     8 - dTSv1
%     9 - dTSv2
%     10 - dTLv1
%     11 - dTLv2
%     12 - dTLv3
%     13 - dTLv4
%     14 - dT La
%   - Global/regional parameters (89)
%     15-21  - VWall  (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4,LA)
%     22-29  - AmRef  (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4,RV,LA)
%     30-37  - SfAct  (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4,RV,LA)
%     38-43  - vMax   (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4)
%     44-51  - TR     (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4,RV,LA)
%     52-57  - TD     (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4)
%     58-64  - Ls0Pas (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4,LA)
%     65-70  - SfPas  (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4)
%     71-77  - k1     (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4,LA)
%     78-83  - dLsPas (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4)
%     84     - LenSE  (global)
%     85-91  - ADO    (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4,LA)
%     92-97  - LDAD   (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4)
%     98-103 - LDCI   (Sv1,Sv2,Lv1,Lv2,Lv3,Lv4)

addpath([pwd,'\CircAdapt-master\matlab\src']);

% Load reference P-struct and create 6-segment LV model
load('PRefHF.mat')
SplitMerge('Sv1',2)
SplitMerge('Lv1',4)
save('PRefHF_6segs.mat','P')

%% Define parameter boundaries
par_bounds = zeros(103,4);

% Global hemodynamics
par_bounds(1,1:2) = [2.1,5.5]/(60*1000);
par_bounds(2,1:2) = P.General.p0*[0.75,1.25];
par_bounds(3,1:2) = [0.6083,1.2907];

% ArtVen
k_Sy_lim = GetFt('ArtVen','k','Sy').*[0.8,1.2];
par_bounds(4,:) = [k_Sy_lim(1),k_Sy_lim(3),k_Sy_lim(2),k_Sy_lim(4)];
p0_Sy_lim = GetFt('ArtVen','p0','Sy').*[0.8,1.2];
par_bounds(5,:) = [p0_Sy_lim(1),p0_Sy_lim(3),p0_Sy_lim(2),p0_Sy_lim(4)];
A0_Sy_lim = GetFt('ArtVen','A0','Sy').*[0.8,1.2];
par_bounds(6,:) = [A0_Sy_lim(1),A0_Sy_lim(3),A0_Sy_lim(2),A0_Sy_lim(4)];

% Patch
% Electrical substrate
par_bounds(7,1:2) = [-0.050,0.130];
par_bounds(8:9,1:2) = repmat([-0.030,0.060],2,1); % dT septum
par_bounds(10:13,1:2) = repmat([0,0.120],4,1); % dT LVfw
par_bounds(14,1:2) = [-0.010,0.050];

% Other parameters
par_bounds([15:20,22:27,44:49,52:57,78:83,85:90],1:2) = repmat([0.64,1.44],36,1);
par_bounds([21,28:29,50:51,77,84,91],1:2) = repmat([0.8,1.2],8,1);
par_bounds(30:35,1:2) = repmat([0,1.44],6,1);
par_bounds(36:37,1:2) = repmat([0.5,1.2],2,1);
par_bounds([38:43,92:103],1:2) = repmat([0.4,1.44],18,1);
par_bounds(58:63,1:2) = repmat([0.9025,1.1025],6,1);
par_bounds(64,1:2) = [0.95,1.05];
par_bounds(65:70,1:2) = repmat([0.64,12],6,1);
par_bounds(71:76,1:2) = repmat([0.64,3.6],6,1);

%% Parameter reference values
ref_par_values = zeros(103,2);

% Global
ref_par_values(1,1) = P.General.q0;
ref_par_values(2,1) = P.General.p0;
ref_par_values(3,1) = P.General.tCycle;

% ArtVen
ref_par_values(4,:) = GetFt('ArtVen','k','Sy')';
ref_par_values(5,:) = GetFt('ArtVen','p0','Sy')';
ref_par_values(6,:) = GetFt('ArtVen','A0','Sy')';

% Patch
% Electrical substrate
ref_par_values(7,1) = P.General.dTauAv;
ref_par_values(8:13,1) = zeros(6,1);
ref_par_values(14,1) = GetFt('Patch','dT','La1');

% VWall
ref_par_values(15:16,1) = repmat(P.Patch.VWall(7),2,1);
ref_par_values(17:20,1) = repmat(P.Patch.VWall(3),4,1);
ref_par_values(21,1) = GetFt('Patch','VWall','La1');

% AmRef
ref_par_values(22:23,1) = repmat(P.Patch.AmRef(7),2,1);
ref_par_values(24:27,1) = repmat(P.Patch.AmRef(3),4,1);
ref_par_values(28,1) = GetFt('Patch','AmRef','Rv1');
ref_par_values(29,1) = GetFt('Patch','AmRef','La1');

% SfAct
ref_par_values(30:35,1) = repmat(P.Patch.SfAct(3),6,1);
ref_par_values(36,1) = GetFt('Patch','SfAct','Rv1');
ref_par_values(37,1) = GetFt('Patch','SfAct','La1');

% vMax
ref_par_values(38:43,1) = repmat(P.Patch.vMax(3),6,1);

% TR
ref_par_values(44:49,1) = repmat(P.Patch.TR(3),6,1);
ref_par_values(50,1) = GetFt('Patch','TR','Rv1');
ref_par_values(51,1) = GetFt('Patch','TR','La1');

% TD
ref_par_values(52:57,1) = repmat(P.Patch.TD(3),6,1);

% Ls0Pas
ref_par_values(58:63,1) = repmat(P.Patch.Ls0Pas(3),6,1);
ref_par_values(64,1) = GetFt('Patch','Ls0Pas','La1');

% SfPas
ref_par_values(65:66,1) = repmat(P.Patch.SfPas(7),2,1);
ref_par_values(67:70,1) = repmat(P.Patch.SfPas(3),4,1);

% k1
ref_par_values(71:76,1) = repmat(P.Patch.k1(3),6,1);
ref_par_values(77,1) = GetFt('Patch','k1','La1');

% dLsPas
ref_par_values(78:83,1) = repmat(P.Patch.dLsPas(3),6,1);

% LenSeriesElement
ref_par_values(84,1) = P.Patch.LenSeriesElement(3);

% ADO
ref_par_values(85:90,1) = repmat(P.Patch.ADO(3),6,1);
ref_par_values(91,1) = GetFt('Patch','ADO','La1');

% LDAD
ref_par_values(92:97,1) = repmat(P.Patch.LDAD(3),6,1);

% LDCI
ref_par_values(98:103,1) = repmat(P.Patch.LDCI(3),6,1);

fields1 = [repmat({'General'},3,1);repmat({'ArtVen'},3,1);'General';repmat({'Patch'},96,1)];
fields2 = ['q0';'p0';'tCycle';'k';'p0';'A0';'dTauAv';repmat({'dT'},7,1);repmat({'VWall'},7,1);repmat({'AmRef'},8,1);repmat({'SfAct'},8,1);repmat({'vMax'},6,1);repmat({'TR'},8,1);repmat({'TD'},6,1);repmat({'Ls0Pas'},7,1);repmat({'SfPas'},6,1);repmat({'k1'},7,1);repmat({'dLsPas'},6,1);'LenSeriesElement';repmat({'ADO'},7,1);repmat({'LDAD'},6,1);repmat({'LDCI'},6,1)];

fields_VWall = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'La1'};
fields_AmRef = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'Rv1';'La1'};
fields_SfAct = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'Rv1';'La1'};
fields_vMax = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4'};
fields_TR = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'Rv1';'La1'};
fields_TD = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4'};
fields_Ls0Pas = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'La1'};
fields_SfPas = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4'};
fields_k1 = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'La1'};
fields_dLsPas = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4'};
fields_LenSeriesElement = 'LV';
fields_ADO = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'La1'};
fields_LDAD = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4'};
fields_LDCI = {'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4'};
fields3 = [repmat({'glob'},3,1);repmat({'Sy'},3,1);'glob';'Sv1';'Sv2';'Lv1';'Lv2';'Lv3';'Lv4';'La1';fields_VWall;fields_AmRef;fields_SfAct;fields_vMax;fields_TR;fields_TD;fields_Ls0Pas;fields_SfPas;fields_k1;fields_dLsPas;fields_LenSeriesElement;fields_ADO;fields_LDAD;fields_LDCI];

%% Initialize simulations
input = sobolset(103,'Skip',1000); % 'Skip' value to be increased by 3000000 for each new set of simulations

n_sims = 3000000; % total number of simulations performed
input_sims = net(input,n_sims);
output_model = cell(n_sims,1);

pars_seg = [8:83,85:103]; % indices of segmental parameters, to be changed when evaluating a different parameter set
pars_glob = 84; % indices of global parameters, to be changed when evaluating a different parameter set

n_substeps = 2; % number of substeps from reference to simulation (for computational stability)

% Translate domain [0,1] to parameter domain
for ind = 1:103
    % Use domain [0.8,1.2] for ArtVen parameters
    if ind>=4 && ind<=6
        input_sims(:,ind) = 0.8+(input_sims(:,ind)*(1.2-0.8));
    else
        input_sims(:,ind) = par_bounds(ind,1)+(input_sims(:,ind)*(par_bounds(ind,2)-par_bounds(ind,1)));
    end
end

% For parameters with a relative parameter value (mf = multiplication
% factor), set absolute parameter value in input_sims
pars_mf = 15:103;
for ind = 1:length(pars_mf)
    input_sims(:,pars_mf(ind)) = input_sims(:,pars_mf(ind)).*ref_par_values(pars_mf(ind),1);
end

% Run all simulations in parallel (calling RunSobolSimulation.m)
parfor n = 1:n_sims
    input_pars = input_sims(n,:)';
    D = RunSobolSimulation(ref_par_values,input_pars,n_substeps,fields1,fields2,fields3,pars_seg,pars_glob);
    output_model{n,1} = D;
end

% To limit the size of the output data file, save as different suboutputs
n_subout = 30;
for k = 1:n_subout
    eval(['output_model',num2str(k),' = output_model(',num2str(((k-1)*n_sims/n_subout)+1),':',num2str(k*n_sims/n_subout),');']);
end

for k = 1:n_subout
    outfolder = ['\Output_Diaphony_sub',num2str(k),'.mat'];
    outfile = ['output_model',num2str(k)];
    eval('save([pwd,outfolder],outfile)');
end

