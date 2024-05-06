% ---CalculateObjectiveFunction_18.m---
%
% Calculate objective function value based on LV cavity volume and strain,
% using an 18-segment LV model
% 
% Koopsen T.
% Last modified: 11/28/2023

function [OFval,OFval_seg] = CalculateObjectiveFunction_18(P,measurement)

% Extract strain and volume data from measurement
S_mea = measurement.strain;
GS_mea = mean(S_mea,2);
GSmin_mea = min(GS_mea); % GLS
LV_EDV_mea = measurement.LV.EDV;
LV_EF_mea = ((measurement.LV.EDV-measurement.LV.ESV)/measurement.LV.EDV)*100;
time_mea = measurement.time';
temp_res_mea = time_mea(2)-time_mea(1);

% Determine time point of 10% global re-lengthening + 50 ms
% (ind_10rel_50ms)

ind_10rel = find(GS_mea==GSmin_mea,1); % start at GLS
while GS_mea(ind_10rel) < GSmin_mea-(0.1*GSmin_mea)
    ind_10rel = ind_10rel+1;
end

inds_50ms = round(0.050/temp_res_mea);
ind_10rel_50ms = min([length(time_mea),ind_10rel+inds_50ms]);

% Extract strain and volume data from P-struct (P)
% Calculate strain
S_sim = zeros(length(P.t),18);

% Determine index of mitral valve closure (n_MVC)
n_start = find(P.Patch.C(:,1)>=0.05*max(P.Patch.C(:,1)),1);
n_end = n_start + find(P.Patch.C(n_start:end,1)<0.05*max(P.Patch.C(:,1)),1)-2;

ind_max_MV_q = find(P.Valve.q(:,5) == max(P.Valve.q(n_start:n_end,5)),1);
n_MVC = ind_max_MV_q + find(P.Valve.q(ind_max_MV_q:end,5)<0.05*P.Valve.q(ind_max_MV_q,5),1)-1;

for seg = 1:18
    S_seg = (((P.Patch.Ls(:,seg+2)./P.Patch.Ls(n_MVC,seg+2))-1)*100);
    S_sim(:,seg) = [S_seg(n_MVC:end);S_seg(1:n_MVC-1)];
end

GSmin_sim = min(mean(S_sim,2));
S_sim = S_sim.*(GSmin_mea/GSmin_sim); % normalization
LV_EDV_sim = max(P.Cavity.V(:,7))*10^6;
LV_EF_sim = ((max(P.Cavity.V(:,7))-min(P.Cavity.V(:,7)))/max(P.Cavity.V(:,7)))*100;
time_sim = 0:0.002:(size(S_sim,1)-1)*0.002;

% Use earlier cutoff point if ind_10rel_50ms exceeds simulated strain
% duration
if time_mea(ind_10rel_50ms)>time_sim(end)
    ind_10rel_50ms = length(time_mea(time_mea<=time_sim(end)));
end

% Resample simulation to sampling frequency of measurement
for seg = 1:18
    S_sim_rs(:,seg) = interp1(time_sim,S_sim(:,seg),time_mea(1:ind_10rel_50ms));
end

% Calculate strain rate of measurement
SR_mea = zeros(size(S_mea,1)-1,18);
for seg = 1:18
    SR_mea(:,seg) = diff(S_mea(:,seg))./temp_res_mea;
end

% Calculate strain rate of simulation
SR_sim = zeros(size(S_sim_rs,1)-1,18);
for seg = 1:18
    SR_sim(:,seg) = diff(S_sim_rs(:,seg))./temp_res_mea;
end

% Calculate objective function value
n_points = ind_10rel_50ms;

% Uncertainty: Otterstad et al. (1997)
OF_weight_LV_EDV = 0.13*LV_EDV_mea;
OF_weight_LV_EF = 0.14*LV_EF_mea;

% Volume errors
OF_LV_EDV = ((LV_EDV_sim-LV_EDV_mea)/OF_weight_LV_EDV)^2;
OF_LV_EF = ((LV_EF_sim-LV_EF_mea)/OF_weight_LV_EF)^2;

% If mean left atrial pressure (mLAP) exceeds 25 mmHg, set infinite OF value
OF_mLAP = 0;
if mean(P.Node.p(:,5)/133.322) > 25
    OF_mLAP = Inf;
end

% Total objective function value (without strain and strain rate, to be added below)
OFval = OF_LV_EDV + OF_LV_EF + OF_mLAP;
% Segmental objective function value
OFval_seg = zeros(1,18);

% Strain and strain rate errors
for seg = 1:18
    OFval_seg(seg) = (1/n_points)*sum((1/2*(S_sim_rs(1:n_points,seg)-S_mea(1:n_points,seg))).^2) + (1/(n_points-1))*sum((1/20*(SR_sim(1:n_points-1,seg)-SR_mea(1:n_points-1,seg))).^2);
end

OFval = OFval + sum(OFval_seg);