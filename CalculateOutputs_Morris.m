% ---CalculateOutputs_Morris.m---
%
% Calculate scalar model outputs of Morris Screening Method
%
% Koopsen T.
% Last modified: 11/27/2023

function [Emin,Eminsys,dEpost,dEpre,dEej,Emeansys,Emeanej,ERminsys,ERminej,ERmeansys,ERmeanej,tsh10,tsh50,tsh90,trel10,trel50] = CalculateOutputs_Morris(P)

% Input arguments:
% P         - CircAdapt P-struct
%
% Output arguments:
% Emin      - peak strain
% Eminsys   - peak systolic strain
% dEpost    - post-systolic shortening
% dEpre     - pre-ejection stretch
% dEej      - ejection stretch
% Emeansys  - mean systolic strain
% Emeanej   - mean ejection strain
% ERminsys  - peak systolic strain rate
% ERminej   - peak ejection strain rate
% ERmeansys - mean systolic strain rate
% ERmeanej  - mean ejection strain rate
% tsh10     - time to 10% shortening
% tsh50     - time to 50% shortening
% tsh90     - time to 90% shortening
% trel10    - time to 10% re-lengthening
% trel50    - time to 50% re-lengthening

% Step 1: Calculate strain in 6-segment LV model

cycletime = P.t-P.t(1); % Cycle starting at t=0

[nAVO,nAVC,nMVO,nMVC] = ValveEvents(P,1,'L'); % Determine valve timings

% Calculate strain in LV segments of interest
segs_all = [1,2,3,4,5,6];
E_all = zeros(length(cycletime),6);
inds = [3,4,5,6,7,8];
for segs = 1:length(segs_all)
    fstr = ((P.Patch.Ls(:,inds(segs))./P.Patch.Ls(nMVC,inds(segs)))-1)*100;
    E_all(:,segs_all(segs)) = [fstr(nMVC:end);fstr(1:nMVC-1)];
end

% Define valve timings (indices) with respect to MVC (t=0)
nAVO = nAVO-nMVC+1;
nAVC = nAVC-nMVC+1;
nMVO = nMVO-nMVC+1;
nMVC = 1;

% Step 2: Calculate global strain indices

E_glob = mean(E_all(:,1:6),2);

E_glob_max = max(E_glob);
if length(E_glob_max)>1
    E_glob_max = E_glob_max(1);
end
E_glob_min = min(E_glob);
if length(E_glob_min)>1
    E_glob_min = E_glob_min(1);
end

% Duplicate global strain curve for index calculation
E_glob_2 = [E_glob;E_glob];
ind = find(E_glob_2 == E_glob_min,1);

% 10% global re-lengthening
while E_glob_2(ind)<E_glob_min + 0.10*(E_glob_max-E_glob_min)
    ind = ind+1;
end
% Use time point which is closest to 10% re-lengthening:
if abs(E_glob_2(ind)-(E_glob_min+0.10*(E_glob_max-E_glob_min))) > abs(E_glob_2(ind-1)-(E_glob_min+0.10*(E_glob_max-E_glob_min)))
    ind = ind-1;
end
ind_10rel_glob = ind;

% 50% global re-lengthening
while E_glob_2(ind)<E_glob_min + 0.50*(E_glob_max-E_glob_min)
    ind = ind+1;
end
% Use time point which is closest to 50% re-lengthening
if abs(E_glob_2(ind)-(E_glob_min+0.50*(E_glob_max-E_glob_min))) > abs(E_glob_2(ind-1)-(E_glob_min+0.50*(E_glob_max-E_glob_min)))
    ind = ind-1;
end
ind_50rel_glob = ind;

if ind_50rel_glob>length(E_glob)
    ind_50rel_glob = ind_50rel_glob-length(E_glob);
end

% Calculate total shortening for all segments (until 50% global re-lengthening)
tot_sh = zeros(1,6);
for seg = 1:6
    E_seg = E_all(nMVC:ind_50rel_glob,seg);

    for m = 1:length(E_seg)-1
        if E_seg(m+1)<E_seg(m)
            tot_sh(seg) = tot_sh(seg)+(E_seg(m+1)-E_seg(m));
        end
    end
end

% Step 3: Calculate strain and strain rate indices

nCol      = size(E_all,2)-1;
Emin      = zeros(1,nCol);
Eminsys   = zeros(1,nCol);
dEpost    = zeros(1,nCol);
dEpre     = zeros(1,nCol);
dEej      = zeros(1,nCol);
Emeansys  = zeros(1,nCol);
Emeanej   = zeros(1,nCol);
ERminsys  = zeros(1,nCol);
ERminej   = zeros(1,nCol);
ERmeansys = zeros(1,nCol);
ERmeanej  = zeros(1,nCol);
tsh10     = zeros(1,nCol);
tsh50     = zeros(1,nCol);
tsh90     = zeros(1,nCol);
trel10    = zeros(1,nCol);
trel50    = zeros(1,nCol);

% Strain rate
ER_all = zeros(size(E_all,1)-1,nCol);
for k = 1:nCol
    ER_all(:,k) = diff(E_all(:,k));
end

% Peak strain (Emin) and index of peak strain (indEmin)
indEmin = zeros(1,6);
for i = 1:nCol
    E_p = E_all(:,i);
    Emin(i) = min(E_p);
    indEmin(i) = find(E_p == min(E_p),1);
end

% Total re-lengthening after reaching peak strain
tot_rel = zeros(1,6);
for seg = 1:6
    E_p = E_all(:,seg);
    for m = indEmin(seg):length(E_p)-1
        if E_p(m+1)>E_p(m)
            tot_rel(seg) = tot_rel(seg) + (E_p(m+1)-E_p(m));
        end
    end
    if E_p(1)>E_p(end)
        tot_rel(seg) = tot_rel(seg) + (E_p(1)-E_p(end));
    end
end

% Peak systolic strain (Eminsys)
for i = 1:nCol
    E_p = E_all(:,i);
    Eminsys(i) = min(E_p(nMVC:nAVC));
end

% Post-systolic shortening (dEpost)
for i = 1:nCol
    E_p = E_all(:,i);
    dEpost(i) = 0;
    if ind_50rel_glob < nAVC
        E_p_2 = [E_p;E_p];
        for m = nAVC:(ind_50rel_glob+length(E_p)-1)
            if E_p_2(m+1)<E_p_2(m)
                dEpost(i) = dEpost(i)+(E_p_2(m+1)-E_p_2(m));
            end
        end
    else
        for m = nAVC:ind_50rel_glob-1
            if E_p(m+1)<E_p(m)
                dEpost(i) = dEpost(i)+(E_p(m+1)-E_p(m));
            end
        end
    end
end

% Pre-ejection stretch (dEpre)
for i = 1:nCol
    E_p = E_all(:,i);
    for m = nMVC:nAVO-1
        if E_p(m+1)>E_p(m)
            dEpre(i) = dEpre(i)+(E_p(m+1)-E_p(m));
        end
    end
end

% Ejection stretch (dEej)
for i = 1:nCol
    E_p = E_all(:,i);
    for m = nAVO:nAVC-1
        if E_p(m+1)>E_p(m)
            dEej(i) = dEej(i)+(E_p(m+1)-E_p(m));
        end
    end
end

% Mean systolic strain (Emeansys)
for i = 1:nCol
    E_p = E_all(:,i);
    Emeansys(i) = mean(E_p(nMVC:nAVC));
end

% Mean ejection strain (Emeanej)
for i = 1:nCol
    E_p = E_all(:,i);
    Emeanej(i) = mean(E_p(nAVO:nAVC));
end

% Peak systolic strain rate (ERminsys)
for i = 1:nCol
    ERp = ER_all(:,i);
    ERminsys(i) = min(ERp(nMVC:nAVC));
end

% Peak ejection strain rate (ERminej)
for i = 1:nCol
    ERp = ER_all(:,i);
    ERminej(i) = min(ERp(nAVO:nAVC));
end

% Mean systolic strain rate (ERmeansys)
for i = 1:nCol
    ERp = ER_all(:,i);
    ERmeansys(i) = mean(ERp(nMVC:nAVC));
end

% Mean ejection strain rate (ERmeanej)
for i = 1:nCol
    ERp = ER_all(:,i);
    ERmeanej(i) = mean(ERp(nAVO:nAVC));
end

% Time to 10% shortening (tsh10)
for i = 1:nCol
    E_p = E_all(:,i);
    sh_seg = 0;
    ind = nMVC;
    while sh_seg > 0.1*tot_sh(i)
        if E_p(ind+1)<E_p(ind)
            sh_seg = sh_seg + (E_p(ind+1)-E_p(ind));
        end
        ind = ind+1;
    end
    diff_sh1 = abs(sh_seg-0.1*tot_sh(i));
    diff_sh2 = abs((sh_seg+(E_p(ind-1)-E_p(ind)))-0.1*tot_sh(i));
    
    t1 = (ind-1-nMVC)*P.General.Dt;
    tsh10(i) = t1 + (diff_sh2/(diff_sh1+diff_sh2))*P.General.Dt;
end

% Time to 50% shortening (tsh50)
for i = 1:nCol
    E_p = E_all(:,i);
    sh_seg = 0;
    ind = nMVC;
    while sh_seg > 0.5*tot_sh(i)
        if E_p(ind+1)<E_p(ind)
            sh_seg = sh_seg + (E_p(ind+1)-E_p(ind));
        end
        ind = ind+1;
    end
    diff_sh1 = abs(sh_seg-0.5*tot_sh(i));
    diff_sh2 = abs((sh_seg+(E_p(ind-1)-E_p(ind)))-0.5*tot_sh(i));
    
    t1 = (ind-1-nMVC)*P.General.Dt;
    tsh50(i) = t1 + (diff_sh2/(diff_sh1+diff_sh2))*P.General.Dt;
end

% Time to 90% shortening (tsh90)
for i = 1:nCol
    E_p = E_all(:,i);
    sh_seg = 0;
    ind = nMVC;
    while sh_seg > 0.9*tot_sh(i)
        if E_p(ind+1)<E_p(ind)
            sh_seg = sh_seg + (E_p(ind+1)-E_p(ind));
        end
        ind = ind+1;
    end
    diff_sh1 = abs(sh_seg-0.9*tot_sh(i));
    diff_sh2 = abs((sh_seg+(E_p(ind-1)-E_p(ind)))-0.9*tot_sh(i));
    
    t1 = (ind-1-nMVC)*P.General.Dt;
    tsh90(i) = t1 + (diff_sh2/(diff_sh1+diff_sh2))*P.General.Dt;
end

% Time to 10% re-lengthening (trel10)
for i = 1:nCol
    E_p = E_all(:,i);
    st_seg = 0;
    ind = indEmin(i);
    while st_seg < 0.1*tot_rel(i)
        if E_p(ind+1)>E_p(ind)
            st_seg = st_seg + (E_p(ind+1)-E_p(ind));
        end
        ind = ind+1;
    end
    diff_st1 = abs(st_seg-0.1*tot_rel(i));
    diff_st2 = abs((st_seg-(E_p(ind)-E_p(ind-1)))-0.1*tot_rel(i));

    t1 = (ind-1-nMVC)*P.General.Dt;
    trel10(i) = t1 + (diff_st2/(diff_st1+diff_st2))*P.General.Dt;
end

% Time to 50% re-lengthening (trel50)
for i = 1:nCol
    E_p = E_all(:,i);
    st_seg = 0;
    ind = indEmin(i);
    while st_seg < 0.5*tot_rel(i)
        if E_p(ind+1)>E_p(ind)
            st_seg = st_seg + (E_p(ind+1)-E_p(ind));
        end
        ind = ind+1;
    end
    diff_st1 = abs(st_seg-0.5*tot_rel(i));
    diff_st2 = abs((st_seg-(E_p(ind)-E_p(ind-1)))-0.5*tot_rel(i));

    t1 = (ind-1-nMVC)*P.General.Dt;
    trel50(i) = t1 + (diff_st2/(diff_st1+diff_st2))*P.General.Dt;
end

