% ---RunSobolSimulation.m---
%
% Evaluate Sobol simulation
%
% Koopsen T.
% Last modified: 11/27/2023

function D = RunSobolSimulation(ref_par_values,sim,n_substeps,fields1,fields2,fields3,pars_seg,pars_glob)
% Input arguments:
% - ref_par_values - reference parameter values
% - sim            - parameter values of simulation to be evaluated
% - n_substeps     - number of substeps from reference to current simulation (for computational stability)
% - fields1        - list of 1st field to set parameter values
% - fields2        - list of 2nd field to set parameter values
% - fields3        - list of 3rd field to set parameter values
% - pars_seg       - segmental parameter indices
% - pars_glob      - global parameter indices

% Output arguments:
% - D              - data struct with all relevant output

load('PRefHF_6segs.mat','P')

% Set all parameter values
for substep = 1:n_substeps
    if substep>1
        % Check stability
        isStable = CA.getIsStable();
        if isStable ~= 1
            break
        end
    end
    
    sim_point = ref_par_values(:,1)+(substep*((sim-ref_par_values(:,1))/n_substeps));
    sim_val4 = (1+(substep*(sim(4)-1)/n_substeps))*ref_par_values(4,:)';
    sim_val5 = (1+(substep*(sim(5)-1)/n_substeps))*ref_par_values(5,:)';
    sim_val6 = (1+(substep*(sim(6)-1)/n_substeps))*ref_par_values(6,:)';

    % Global hemodynamics
    for par = 1:3
        eval(['P.General.',fields2{par},'=sim_point(par);'])
    end
    % ArtVen parameters
    PutFt(fields1{4},fields2{4},fields3{4},sim_val4)
    PutFt(fields1{5},fields2{5},fields3{5},sim_val5)
    PutFt(fields1{6},fields2{6},fields3{6},sim_val6)
    % dTauAv
    P.General.dTauAv = sim_point(7);
    % Local patch parameters
    for lp = 1:length(pars_seg)
        PutFt(fields1{pars_seg(lp)},fields2{pars_seg(lp)},fields3{pars_seg(lp)},sim_point(pars_seg(lp)))
    end
    % Global patch parameters
    for gp = 1:length(pars_glob)
        if pars_glob(gp) == 15 % VWall
            PutFt(fields1{pars_glob(gp)},fields2{pars_glob(gp)},{'Lv1','Lv2','Lv3','Lv4'},1.4836*sim_point(pars_glob(gp)))
            PutFt(fields1{pars_glob(gp)},fields2{pars_glob(gp)},{'Sv1','Sv2'},sim_point(pars_glob(gp)))
        elseif pars_glob(gp) == 65 % SfPas
            PutFt(fields1{pars_glob(gp)},fields2{pars_glob(gp)},{'Lv1','Lv2','Lv3','Lv4'},1.0102*sim_point(pars_glob(gp)))
            PutFt(fields1{pars_glob(gp)},fields2{pars_glob(gp)},{'Sv1','Sv2'},sim_point(pars_glob(gp)))
        else
            PutFt(fields1{pars_glob(gp)},fields2{pars_glob(gp)},{'Sv1','Sv2','Lv1','Lv2','Lv3','Lv4'},sim_point(pars_glob(gp)))
        end
    end
    
    RunCPP
end

% Postprocessing: calculate strain, volumes and pressures
% (1) 6-segment strain (D.strain)
% (2) LV volumes (EDV,ESV) (D.LV)
% (3) RV volumes (EDV,ESV) (D.RV)
% (4) LA volumes (EDV,ESV) (D.LA)
% (5) Mean LA pressure (D.mLAP)
% (6) Valve timings (AVO,AVC,MVO,MVC) (D.valves)

% Check stability
isStable = CA.getIsStable();
if isStable == 1
    try
        % Fill D-struct
        % (7) Valve timings (MVC) (D.valves)
        % MVC
        n_start = find(P.Patch.C(:,1)>=0.05*max(P.Patch.C(:,1)),1);
        n_end = n_start + find(P.Patch.C(n_start:end,1)<0.05*max(P.Patch.C(:,1)),1)-2;
               
        ind_max_MV_q = find(P.Valve.q(:,5) == max(P.Valve.q(n_start:n_end,5)),1);
        n_MVC = ind_max_MV_q + find(P.Valve.q(ind_max_MV_q:end,5)<0.05*P.Valve.q(ind_max_MV_q,5),1)-1;
        
        D.valves.MVC = n_MVC;

        % (2) 6-segment LV strain (D.strain)
        strain_all = zeros(length(P.t),6);
        for segs = 1:6
            cstr = ((P.Patch.Ls(:,segs+2)./P.Patch.Ls(n_MVC,segs+2))-1)*100;
            strain_all(:,segs) = [cstr(n_MVC:end);cstr(1:n_MVC-1)];
        end
        D.strain = strain_all;

        % (3) LV volumes (D.LV)
        D.LV.EDV = max(P.Cavity.V(:,7))*10^6;
        D.LV.ESV = min(P.Cavity.V(:,7))*10^6;

        % (4) RV volumes (D.RV)
        D.RV.EDV = max(P.Cavity.V(:,8))*10^6;
        D.RV.ESV = min(P.Cavity.V(:,8))*10^6;

        % (5) LA volumes (D.LA)
        %D.LA.EDV = max(P.Cavity.V(:,5))*10^6;
        %D.LA.ESV = min(P.Cavity.V(:,5))*10^6;

        % (6) Mean LA pressure (D.mLAP)
        D.mLAP = mean(P.Node.p(:,5));
    catch
        %disp('Strain calculation failed')
        D = NaN;
    end
else
    %disp('Simulation not stable')
    D = NaN;
end

