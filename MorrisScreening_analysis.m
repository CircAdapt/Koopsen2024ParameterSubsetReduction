% ---MorrisScreening_analysis.m---
%
% Evaluates trajectories of Morris Screening Method
%
% Koopsen T.
% Last modified: 11/27/2023

function [delta_strain_indices_LV,delta_LVEDV,delta_LVSV] = MorrisScreening_analysis(P_start,all_lims,fields1,fields2,fields3,nsubsteps_sp,nsubsteps_tr,starting_point_0,traj,savefolder,n)
% Input arguments:
% P_start                      - P-struct used as starting point for all simulations
% all_lims                     - all parameter limits (defined in Morris_screening_input.m)
% fields1                      - list of 1st field to set parameter values
% fields2                      - list of 2nd field to set parameter values
% fields3                      - list of 3rd field to set parameter values
% nsubsteps_sp                 - number of substeps from P_start to starting point of trajectory
% nsubsteps_tr                 - number of substeps between different points within trajectory
% starting_point_0             - normalized starting points of trajectories
% traj                         - trajectory definition (includes sequence of parameter changes)
% savefolder                   - folder in which output is saved
% n                            - trajectory number
%
% Output arguments:
% delta_strain_indices_LV      - absolute change of LV strain indices
% delta_LVEDV                  - absolute change of LV end-diastolic volume
% delta_LVSV                   - absolute change of LV stroke volume

% Load P_start
load(P_start,'P')

% Create output matrices
delta_strain_indices_LV = zeros(174,6,16);

delta_LVEDV = zeros(174,1);
delta_LVSV = zeros(174,1);

strain_indices_names = {'Emin','Eminsys','dEpost','dEpre','dEej','Emeansys','Emeanej','ERminsys','ERminej','ERmeansys','ERmeanej','tsh10','tsh50','tsh90','trel10','trel50'};
for k = 1:length(strain_indices_names)
    eval(['delta_',strain_indices_names{k},' = zeros(174,6);'])
end

% Subdivide LV wall into 6 patches (6-segment model)
% Names used in this script:                              Sv1 - AntSept, Sv2 - Sept, Lv1 - Ant, Lv2 - Inf, Lv3 - Lat, Lv4 - Post
% Names used in Koopsen et al. (2024) Biomed. Eng. Onl.:  S1  - AntSept, S2  - Sept, LV1 - Ant, LV2 - Inf, LV3 - Lat, LV4 - Post
SplitMerge('Sv1',2);
SplitMerge('Lv1',4);

% Define direction of electrical activation (LBBB-like activation)
Early = 'Sv1';
Late = 'Lv4';
Medium_Early1 = 'Lv1';
Medium_Late1 = 'Lv3';
Medium_Early2 = 'Sv2';
Medium_Late2 = 'Lv2';

% Run CircAdapt model
RunCPP

% Check stability of simulation
isStable = CA.getIsStable();

for substep = 1:nsubsteps_sp
    if isStable == 0 % if simulation is not stable, exit
        break
    end
    
    % Define intermediate starting point based on the number of substeps taken
    starting_point = (substep*((starting_point_0-(3.5/7))/nsubsteps_sp))+(3.5/7);

    for par = 1:174
        % Not all parameters have a symmetrical upper and lower limit compared to their reference value,
        % therefore, define starting points specifically for these parameters
        if par==3 % tCycle
            starting_point(par) = (substep*((starting_point_0(par)-0.3294)/nsubsteps_sp))+0.3294;
        elseif par==27 % dTauAv
            starting_point(par) = (substep*((starting_point_0(par)-0.3333)/nsubsteps_sp))+0.3333;
        elseif par==29 % tau2
            starting_point(par) = (substep*((starting_point_0(par)-0)/nsubsteps_sp))+0;
        elseif par==34 % dT La
            starting_point(par) = (substep*((starting_point_0(par)-0.1667)/nsubsteps_sp))+0.1667;
        elseif par==55 || par==62 || par==63 || par==64 || par==66 || par==68 || par==70 || par==156 || par==158 || par==160 || par==166 || par==168 || par==170
            % SfAct GO, SfAct RV, SfAct LA, SfAct RA, vMax I (3), LDAD I (3), LDCI I (3)
            starting_point(par) = (substep*((starting_point_0(par)-0.7143)/nsubsteps_sp))+0.7143;
        elseif par==56 || par==58 || par==60 % SfAct I (3)
            starting_point(par) = (substep*((starting_point_0(par)-0.8333)/nsubsteps_sp))+0.8333;
        elseif par==106 || par==108 || par==110 % SfPas I (3)
            starting_point(par) = (substep*((starting_point_0(par)-0.0217)/nsubsteps_sp))+0.0217;
        elseif par==116 || par==118 || par==120 % k1 I (3)
            starting_point(par) = (substep*((starting_point_0(par)-0.0909)/nsubsteps_sp))+0.0909;
        end
    end

    for par = 1:174
        isStable = CA.getIsStable(); % check stability of simulation
        if isStable == 1
            % Set parameter to target value
            if par<=3
                eval(['P.General.',fields2{par},'=all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1));'])
            elseif par>=4 && par<=13
                PutFt(fields1{par},fields2{par},fields3{par},[all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1));all_lims(par,3)+starting_point(par)*(all_lims(par,4)-all_lims(par,3))])
            elseif par>13 && par<=26
                PutFt(fields1{par},fields2{par},fields3{par},all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))
            elseif par==27
                eval(['P.General.',fields2{par},'=all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1));'])
            elseif par==28
                Onset_sep = max([0,all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))]);
                Onset_RV = -1*min([0,all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))]);
                PutFt('Patch','dT',Early,Onset_sep)
                PutFt('Patch','dT','Rv1',Onset_RV)
            elseif par==29
                Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                PutFt('Patch','dT',Late,Onset_sep+(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))))
            elseif par==30
                Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                PutFt('Patch','dT',Medium_Early1,Onset_sep+((all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
            elseif par==31
                Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                PutFt('Patch','dT',Medium_Late1,Onset_sep+((all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
            elseif par==32
                Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                PutFt('Patch','dT',Medium_Early2,Onset_sep+((all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
            elseif par==33
                Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                PutFt('Patch','dT',Medium_Late2,Onset_sep+((all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
            elseif par==34
                PutFt('Patch','dT','La1',all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))
            elseif par==35 || par==45 || par==55 || par==65 || par==75 || par==85 || par==95 || par==105 || par==115 || par==125 || par==135 || par==145 || par==155 || par==165
                GO = all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1));
                if par==35 || par==45
                    % septum
                    for seg = 1:2
                        PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},0.12603044*GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                    end
                    % LV free wall
                    for seg = 3:6
                        PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},0.18698478*GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                    end
                elseif par==105
                    % septum
                    for seg = 1:2
                        PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},22086*GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                    end
                    % LV free wall
                    for seg = 3:6
                        PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},22311*GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                    end
                else
                    for seg = 1:6
                        PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                    end
                end
            elseif mod(par,10)==6 || mod(par,10)==7 || mod(par,10)==8 || mod(par,10)==9 || mod(par,10)==0 || mod(par,10)==1
                parstr = num2str(par-6);
                GO_par = str2num([parstr(1:length(parstr)-1),'5']);
                GO = all_lims(GO_par,1)+starting_point(GO_par)*(all_lims(GO_par,2)-all_lims(GO_par,1));
                if par>=36 && par<=37
                    PutFt(fields1{par},fields2{par},fields3{par},0.12603044*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                elseif par>=38 && par<=41
                    PutFt(fields1{par},fields2{par},fields3{par},0.18698478*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                elseif par>=46 && par<=47
                    PutFt(fields1{par},fields2{par},fields3{par},0.12603044*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                elseif par>=48 && par<=51
                    PutFt(fields1{par},fields2{par},fields3{par},0.18698478*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                elseif par>=106 && par<=107
                    PutFt(fields1{par},fields2{par},fields3{par},22086*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                elseif par>=108 && par<=111
                    PutFt(fields1{par},fields2{par},fields3{par},22311*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                else
                    PutFt(fields1{par},fields2{par},fields3{par},GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                end
            else
                PutFt(fields1{par},fields2{par},fields3{par},all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)));
            end
            RunCPP
        else
            break
        end
    end
end

isStable = CA.getIsStable();
if isStable == 1
    % Try to calculate all LV strain and volume indices
    try
        [Emin_1,Eminsys_1,dEpost_1,dEpre_1,dEej_1,Emeansys_1,Emeanej_1,ERminsys_1,ERminej_1,ERmeansys_1,ERmeanej_1,tsh10_1,tsh50_1,tsh90_1,trel10_1,trel50_1] = CalculateOutputs_Morris(P);
        LVEDV_1 = max(P.Cavity.V(:,7))*10^6;
        LVSV_1 = (max(P.Cavity.V(:,7))-min(P.Cavity.V(:,7)))*10^6;
    catch
        isStable = 0;
    end
end

if isStable == 1
    finished = 0;

    for step = 1:174
        if finished+1 == step
            par = traj(step);
            for substep = 1:nsubsteps_tr
                if starting_point(par)<4/7
                    starting_point(par) = starting_point(par)+((4/7)/nsubsteps_tr);
                else
                    starting_point(par) = starting_point(par)-((4/7)/nsubsteps_tr);
                end
                if par<=3
                    eval(['P.General.',fields2{par},'=all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1));'])
                elseif par>=4 && par<=13
                    PutFt(fields1{par},fields2{par},fields3{par},[all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1));all_lims(par,3)+starting_point(par)*(all_lims(par,4)-all_lims(par,3))])
                elseif par>13 && par<=26
                    PutFt(fields1{par},fields2{par},fields3{par},all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))
                elseif par==27
                    eval(['P.General.',fields2{par},'=all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1));'])
                elseif par==28
                    Onset_sep = max([0,all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))]);
                    Onset_RV = -1*min([0,all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))]);
                    PutFt('Patch','dT',Early,Onset_sep)
                    PutFt('Patch','dT',Late,Onset_sep+(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1))))
                    PutFt('Patch','dT',Medium_Early1,Onset_sep+((all_lims(30,1)+starting_point(30)*(all_lims(30,2)-all_lims(30,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
                    PutFt('Patch','dT',Medium_Late1,Onset_sep+((all_lims(31,1)+starting_point(31)*(all_lims(31,2)-all_lims(31,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
                    PutFt('Patch','dT',Medium_Early2,Onset_sep+((all_lims(32,1)+starting_point(32)*(all_lims(32,2)-all_lims(32,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
                    PutFt('Patch','dT',Medium_Late2,Onset_sep+((all_lims(33,1)+starting_point(33)*(all_lims(33,2)-all_lims(33,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
                    PutFt('Patch','dT','Rv1',Onset_RV)
                elseif par==29
                    Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                    PutFt('Patch','dT',Late,Onset_sep+(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))))
                    PutFt('Patch','dT',Medium_Early1,Onset_sep+((all_lims(30,1)+starting_point(30)*(all_lims(30,2)-all_lims(30,1)))*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))))
                    PutFt('Patch','dT',Medium_Late1,Onset_sep+((all_lims(31,1)+starting_point(31)*(all_lims(31,2)-all_lims(31,1)))*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))))
                    PutFt('Patch','dT',Medium_Early2,Onset_sep+((all_lims(32,1)+starting_point(32)*(all_lims(32,2)-all_lims(32,1)))*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))))
                    PutFt('Patch','dT',Medium_Late2,Onset_sep+((all_lims(33,1)+starting_point(33)*(all_lims(33,2)-all_lims(33,1)))*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))))
                elseif par==30
                    Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                    PutFt('Patch','dT',Medium_Early1,Onset_sep+((all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
                elseif par==31
                    Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                    PutFt('Patch','dT',Medium_Late1,Onset_sep+((all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
                elseif par==32
                    Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                    PutFt('Patch','dT',Medium_Early2,Onset_sep+((all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
                elseif par==33
                    Onset_sep = max([0,all_lims(28,1)+starting_point(28)*(all_lims(28,2)-all_lims(28,1))]);
                    PutFt('Patch','dT',Medium_Late2,Onset_sep+((all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))*(all_lims(29,1)+starting_point(29)*(all_lims(29,2)-all_lims(29,1)))))
                elseif par==34
                    PutFt('Patch','dT','La1',all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)))
                elseif par==35 || par==45 || par==55 || par==65 || par==75 || par==85 || par==95 || par==105 || par==115 || par==125 || par==135 || par==145 || par==155 || par==165
                    GO = all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1));
                    if par==35 || par==45
                        for seg = 1:2
                            PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},0.12603044*GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                        end
                        for seg = 3:6
                            PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},0.18698478*GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                        end
                    elseif par==105
                        for seg = 1:2
                            PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},22086*GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                        end
                        for seg = 3:6
                            PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},22311*GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                        end
                    else
                        for seg = 1:6
                            PutFt(fields1{par+seg},fields2{par+seg},fields3{par+seg},GO*(all_lims(par+seg,1)+starting_point(par+seg)*(all_lims(par+seg,2)-all_lims(par+seg,1))));
                        end
                    end
                elseif mod(par,10)==6 || mod(par,10)==7 || mod(par,10)==8 || mod(par,10)==9 || mod(par,10)==0 || mod(par,10)==1
                    parstr = num2str(par-6);
                    GO_par = str2num([parstr(1:length(parstr)-1),'5']);
                    GO = all_lims(GO_par,1)+starting_point(GO_par)*(all_lims(GO_par,2)-all_lims(GO_par,1));
                    if par>=36 && par<=37
                        PutFt(fields1{par},fields2{par},fields3{par},0.12603044*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                    elseif par>=38 && par<=41
                        PutFt(fields1{par},fields2{par},fields3{par},0.18698478*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                    elseif par>=46 && par<=47
                        PutFt(fields1{par},fields2{par},fields3{par},0.12603044*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                    elseif par>=48 && par<=51
                        PutFt(fields1{par},fields2{par},fields3{par},0.18698478*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                    elseif par>=106 && par<=107
                        PutFt(fields1{par},fields2{par},fields3{par},22086*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                    elseif par>=108 && par<=111
                        PutFt(fields1{par},fields2{par},fields3{par},22311*GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                    else
                        PutFt(fields1{par},fields2{par},fields3{par},GO*(all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1))));
                    end
                else
                    PutFt(fields1{par},fields2{par},fields3{par},all_lims(par,1)+starting_point(par)*(all_lims(par,2)-all_lims(par,1)));
                end
                RunCPP
            end
            isStable = CA.getIsStable();
            if isStable == 1
                try
                    [Emin_2,Eminsys_2,dEpost_2,dEpre_2,dEej_2,Emeansys_2,Emeanej_2,ERminsys_2,ERminej_2,ERmeansys_2,ERmeanej_2,tsh10_2,tsh50_2,tsh90_2,trel10_2,trel50_2] = CalculateOutputs_Morris(P);
                    LVEDV_2 = max(P.Cavity.V(:,7))*10^6;
                    LVSV_2 = (max(P.Cavity.V(:,7))-min(P.Cavity.V(:,7)))*10^6;
                    
                    for k = 1:length(strain_indices_names)
                        eval(['delta_',strain_indices_names{k},'(par,:) = abs(',strain_indices_names{k},'_2-',strain_indices_names{k},'_1);'])
                        eval([strain_indices_names{k},'_1 = ',strain_indices_names{k},'_2;'])
                    end

                    delta_LVEDV(par) = abs(LVEDV_2-LVEDV_1);
                    delta_LVSV(par) = abs(LVSV_2-LVSV_1);

                    LVEDV_1 = LVEDV_2;
                    LVSV_1 = LVSV_2;

                    finished = finished+1;
                catch
                    %disp('Trajectory failed')                       
                end
            else
                %disp('Simulation not stable')
            end
        else
        end
    end
    if finished == 174
        for k = 1:length(strain_indices_names)
            eval(['delta_strain_indices_LV(:,:,k) = delta_',strain_indices_names{k},';'])
        end
                
        %files = dir(savefolder);
        save([savefolder,'\Trajectory',num2str(n),'.mat'],'delta_strain_indices_LV','delta_LVEDV','delta_LVSV')

        disp('Trajectory successful')
    else
    end
end
end