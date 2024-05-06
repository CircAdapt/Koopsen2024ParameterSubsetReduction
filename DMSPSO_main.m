% ---DMSPSO_main.m---
%
% Dynamic multi-swarm particle swarm optimization (DMS-PSO) algorithm to 
% personalize CircAdapt to patient-specific LV cavity volume and strain data
%
% Koopsen T.
% Last modified: 11/27/2023

function [global_opt_params,global_opt_OFval] = DMSPSO_main(datafile,par_set,parallel_id,rep_id)

% Input arguments:
% datafile          - data struct containing patient-specific input data
% par_set           - parameter set used for optimization
% parallel_id       - unique optimization ID for parallel evaluations
% rep_id            - unique repetition ID for repeated evaluations
%
% Output arguments:
% global_opt_params - parameter values of global optimum particle (end result of optimization)
% global_opt_OFval  - objective function value of global optimum particle

% Independent of par_set, the following parameter indexing is used in this script:
% 1       - q0
% 2       - p0
% 3       - k_Sy
% 4       - p0_Sy
% 5       - A0_Sy
% 6       - dTauAv
% 7-24    - dT_Lv (18 patches)
% 25      - dT_La
% 26-43   - VWall_Lv (18 patches)
% 44      - VWall_La
% 45-62   - AmRef_Lv (18 patches)
% 63      - AmRef_Rv
% 64      - AmRef_La
% 65-82   - SfAct_Lv (18 patches)
% 83      - SfAct_Rv
% 84      - SfAct_La
% 85-102  - vMax_Lv (18 patches)
% 103-120 - TR_Lv (18 patches)
% 121     - TR_Rv
% 122     - TR_La
% 123-140 - TD_Lv (18 patches)
% 141-158 - Ls0Pas_Lv (18 patches)
% 159     - Ls0Pas_La
% 160-177 - SfPas_Lv (18 patches)
% 178-195 - k1_Lv (18 patches)
% 196     - k1_La
% 197-214 - dLsPas_Lv (18 patches)
% 215     - LenSE_Lv_G (global)
% 216-233 - ADO_Lv (18 patches)
% 234     - ADO_La
% 235-252 - LDAD_Lv (18 patches)
% 253-270 - LDCI_Lv (18 patches)

% Add path with CircAdapt C++ code
CircAdaptMatlab = [pwd,'\CircAdapt-master\matlab\src'];
addpath(CircAdaptMatlab);

%% Load datafile
% 'datafile' contains struct 'measurement' with fields:
% 'strain' - array with strain data (Ntime*Nseg)
% 'LV'     - contains two subfields:
%   'EDV'    - LV end-diastolic volume (mL)
%   'ESV'    - LV end-systolic volume (mL)
% 'time'   - array with time data (Ntime*1)
% 'ctmean' - average cycle time of three acquisitions (s)

load([datafile,'.mat'],'measurement');

% Perform a translation from the standard AHA segmentation
% (basAnt, basAntSept, basSept, basInf, basPost, basLat, midAnt,
% midAntSept, midSept, midInf, midPost, midLat, apAnt, apAntSept, apSept,
% apInf, apPost, apLat)
% to the CircAdapt patch sequence
% (basAnt, basInf, basPost, basLat, midAnt, midInf, midPost, midLat, apAnt,
% apInf, apPost, apLat, basAntSept, basSept, midAntSept, midSept,
% apAntSept, apSept)

measurement_copy = measurement;
seq_AHA_to_CA = [1,4,5,6,7,10,11,12,13,16,17,18,2,3,8,9,14,15];
for seg = 1:18
    measurement_copy.strain(:,seg) = measurement.strain(:,seq_AHA_to_CA(seg));
end
measurement = measurement_copy;


%% Define parameter boundaries (par_bounds) for Monte Carlo (MC) simulations
par_bounds = setBounds('PRef_AHA18.mat',par_set,measurement);

%% Save P-struct with cycle time set to measurement value
load('PRef_AHA18_CPP.mat','P');
P.General.tCycle = measurement.ctmean;
RunCPP
save(['PRef_CPP_ctset_',num2str(parallel_id),'.mat'],'P');

%% Perform Monte Carlo (MC) simulations
nsims = 1000; % total number of MC simulations

npars = 270; % total number of parameters (not necessarily all used in optimization, this depends on par_set)

OFvals = zeros(nsims,1); % all objective function values
OFvals_seg = zeros(nsims,18); % all objective function values per segment

rand_sims = rand(nsims,npars); % MC simulations to be evaluated
pars = zeros(nsims,npars); % parameter values

% Run MC simulations
for sim = 1:nsims
    %disp(num2str(sim))
    % Translate domain [0,1] to parameter domain
    pars(sim,:) = par_bounds(:,1) + rand_sims(sim,:)'.*(par_bounds(:,2)-par_bounds(:,1));
    
    % ADDITIONAL RULE 1: Simulations without contractile dysfunction should
    % be included more than in random sampling, therefore:
    % - Draw random number of patches between 1 and 18 and assign
    %   relatively normal SfAct (between 96 and 144 kPa) to these patches
    if contains('SfAct_Lv',par_set)
        seg_seq = randperm(18);
        segs_no_cdysf = seg_seq(1:randi(18));
        pars(sim,segs_no_cdysf+64) = 96000+rand(1,length(segs_no_cdysf))*48000;
    end
    
    % ADDITIONAL RULE 2: Global stiffness should not be too high,
    % therefore:
    % - Draw random number of patches between 1 and 18 and assign
    %   relatively normal SfPas (between 80% and 120% of REF) and k1
    %   (between 8 and 12) to these patches
    if contains('SfPas_Lv',par_set) || contains('k1_Lv',par_set)
        seg_seq = randperm(18);
        segs_low_stiff = seg_seq(1:randi(18));
        if contains('SfPas_Lv',par_set)
            pars(sim,segs_low_stiff(segs_low_stiff<=6)+159) = (0.8*22086)+rand(1,length(segs_low_stiff(segs_low_stiff<=6)))*(0.4*22086);
            pars(sim,segs_low_stiff(segs_low_stiff>6)+159) = (0.8*22311)+rand(1,length(segs_low_stiff(segs_low_stiff>6)))*(0.4*22311);
        end
        if contains('k1_Lv',par_set)
            pars(sim,segs_low_stiff+177) = 8+rand(1,length(segs_low_stiff))*4;
        end
    end
    
    % ADDITIONAL RULE 3: Electrical substrate should not always be too
    % extreme, more synchronous simulations should also be possible,
    % therefore:
    % - Draw random number between 1 and 3 as the 'severity index' of
    %   electrical dyssynchrony, where 1 is least and 3 is most severe, and
    %   scale dT values accordingly
    el_severity = randi(3);
    pars(sim,7:24) = pars(sim,7:24)*(el_severity/3);
    
    % Load P-struct and set parameter values
    P = setSim(['PRef_CPP_ctset_',num2str(parallel_id),'.mat'],par_set,pars(sim,:));
    
    % Modified RunCPP script which is faster (not all output calculated,
    % only output used in objective function)
    RunCPP2
    
    % Calculate objective function value
    try
        % Cavity volumes should not be negative
        if min(P.Cavity.V(:,7).*10^6) <= 0 || min(P.Cavity.V(:,8).*10^6) <= 0
            error('Volume criteria not met')
        end
        
        % Cavity pressures should not be negative
        if min(P.Node.p(:,8)./133.322) <= 0 || min(P.Node.p(:,7)./133.322) <= 0
            error('Pressure criteria not met')
        end
        
        % Calculate objective function value
        [OFval,OFval_seg] = CalculateObjectiveFunction_18(P,measurement);
        
        OFvals(sim) = OFval;
        OFvals_seg(sim,:) = OFval_seg;
    catch
        %disp(['Simulation ', num2str(sim),' did not run'])
        OFvals(sim) = Inf;
        OFvals_seg(sim,:) = zeros(1,18)+Inf;
    end
end

%% Initialize dynamic multi-swarm particle swarm optimization (DMS-PSO) 
% First: pick best MC simulations as starting points for particles

nparticles = 60; % number of particles
nswarms = 20; % number of swarms
max_niter = 2000; % maximum number of DMS-PSO iterations
nsteps_np = 1; % number of steps towards new particle position

% Keep track of best particle objective function value (ever encountered during DMS-PSO) and best
% corresponding parameters
particle_best_OFvals = zeros(1,nparticles);
particle_cur_OFvals = zeros(1,nparticles);
particle_best_params = zeros(npars,nparticles);

% n_GLOB_OF and n_SEGM_OF allow the best particles to be selected based on
% the overall objective function value (GLOB) or the segmental objective function
% value (SEGM). In the current version of the algorithm, as also used in
% Koopsen et al. (2024) Biomed. Eng. Onl. particles are chosen based on the
% best global objective function value, therefore n_GLOB_OF equals the
% total number of particles
n_GLOB_OF = nparticles; % number of particles selected as best overall objective function value (all 18 segments) - GLOB_OF
n_SEGM_OF = nparticles-n_GLOB_OF; % number of particles selected from segmental objective function value (individual segment) - SEGM_OF

% GLOB_OF
best_GLOB_OFinds = zeros(nparticles,1); % nparticles and not n_GLOB_OF since also used for combined segmental simulations
best_GLOB_OFvals = zeros(nparticles,1);

% Find nparticles best overall OF indices and values (not only used for best GLOB_OF but also for combined segmental simulations)
for best = 1:nparticles
    ind = find(OFvals == min(OFvals),1);
    best_GLOB_OFinds(best) = ind;
    best_GLOB_OFvals(best) = OFvals(ind);
    OFvals(ind) = Inf;
end

% SEGM_OF
% best_SEGM_OF_inds = zeros(n_SEGM_OF,18);

% for seg = 1:18
%     for best = 1:n_SEGM_OF
%         ind = find(errs_seg(:,seg) == min(errs_seg(:,seg)),1);
%         best_SEGM_OF_inds(best,seg) = ind;
%         errs_seg(ind,seg) = Inf;
%     end
% end

% To obtain random assignment to swarms, randomize sequences

% GLOB_OF
rand_seq = randperm(n_GLOB_OF);
best_GLOB_OFinds = best_GLOB_OFinds(rand_seq);

% SEGM_OF
% for seg = 1:18
%     rand_seq = randperm(n_SEGM_OF);
%     best_SEGM_OF_inds(:,seg) = best_SEGM_OF_inds(rand_seq,seg);
% end
% rand_seq = randperm(n_SEGM_OF)+n_GLOB_OF; % select 'next best n_SEGM_OF' overall OF simulations
% best_next_GLOB_OF_inds = best_GLOB_OF_inds(rand_seq);

% Run particles and save to workspace

glob_seg_arr = [ones(1,n_GLOB_OF),ones(1,n_SEGM_OF)+1];
rand_seq = randperm(nparticles);
glob_seg_arr = glob_seg_arr(rand_seq);

glob_counter = 0;
seg_counter = 0;

for particle = 1:nparticles
    
    glob_seg = glob_seg_arr(particle);
    if glob_seg == 1
        glob_counter = glob_counter+1;
        par_vec = pars(best_GLOB_OFinds(glob_counter),:);
    elseif glob_seg == 2
        seg_counter = seg_counter+1;
        par_vec = zeros(1,npars);
        % Parameters other than LV patch parameters from best next n_SEGM_OF global simulations
        par_vec([1:6,25,44,63:64,83:84,121:122,159,196,215,234]) = pars(best_next_GLOB_OF_inds(seg_counter),[1:6,25,44,63:64,83:84,121:122,159,196,215,234]);
        for seg = 1:18
            par_vec([seg+6,seg+25,seg+44,seg+84,seg+102,seg+122,seg+140,seg+159,seg+177,seg+196,seg+215,seg+234,seg+252]) = pars(best_SEGM_OF_inds(seg_counter,seg),[seg+6,seg+25,seg+44,seg+84,seg+102,seg+122,seg+140,seg+159,seg+177,seg+196,seg+215,seg+234,seg+252]);
        end
    end
    
    % Set parameter values
    P = setSim(['PRef_CPP_ctset_',num2str(parallel_id),'.mat'],par_set,par_vec);
    
    RunCPP2
    
    if glob_seg == 1
        particle_best_OFvals(particle) = best_GLOB_OFvals(glob_counter);
        particle_cur_OFvals(particle) = best_GLOB_OFvals(glob_counter);
    elseif glob_seg == 2
        try
            if min(P.Cavity.V(:,7).*10^6) <= 0 || min(P.Cavity.V(:,8).*10^6) <= 0
                error('Volume criteria not met')
            end
            
            if min(P.Node.p(:,8)./133.322) <= 0 || min(P.Node.p(:,7)./133.322) <= 0
                error('Pressure criteria not met')
            end
            
            % Calculate objective function
            [OFval,~] = CalculateObjectiveFunction_18(P,measurement);
        
            particle_best_OFvals(particle) = OFval;
            particle_cur_OFvals(particle) = OFval;
        catch
            disp('Particle from segmental errors did not run, trying other combination')
            particle_best_OFvals(particle) = Inf;
            while particle_best_OFvals(particle) == Inf
                for seg = 1:18
                    par_vec([seg+6,seg+25,seg+44,seg+84,seg+102,seg+122,seg+140,seg+159,seg+177,seg+196,seg+215,seg+234,seg+252]) = pars(best_SEGM_OF_inds(randi(n_SEGM_OF),seg),[seg+6,seg+25,seg+44,seg+84,seg+102,seg+122,seg+140,seg+159,seg+177,seg+196,seg+215,seg+234,seg+252]);
                end
                
                % Set parameter values
                P = setSim(['PRef_CPP_ctset_',num2str(parallel_id),'.mat'],par_set,par_vec);
                RunCPP2
                
                try
                    if min(P.Cavity.V(:,7).*10^6) <= 0 || min(P.Cavity.V(:,8).*10^6) <= 0
                        error('Volume criteria not met')
                    end
                    
                    if min(P.Node.p(:,8)./133.322) <= 0 || min(P.Node.p(:,7)./133.322) <= 0
                        error('Pressure criteria not met')
                    end
                    
                    % Calculate objective function
                    [OFval_i,~] = CalculateObjectiveFunction_18(P,measurement);

                    particle_best_OFvals(particle) = OFval_i;
                    particle_cur_OFvals(particle) = OFval_i;
                catch
                    particle_best_OFvals(particle) = Inf;
                    particle_cur_OFvals(particle) = Inf;
                end
            end
        end
    end
    particle_best_params(:,particle) = par_vec';
    save(['P_PSO',num2str(parallel_id),'_Particle',num2str(particle),'.mat'],'P')
end

% Keep track of current particle parameters
particle_cur_params = particle_best_params;

% Define subswarms
swarm_groups = zeros(nparticles/nswarms,nswarms);
for sw = 1:nswarms
    swarm_groups(:,sw) = sw*(nparticles/nswarms)-((nparticles/nswarms)-1):sw*(nparticles/nswarms);
end

% Initialize output
best_swarm_OFval = zeros(1,nswarms); % best objective function value per swarm (used as the swarm's optimum)
best_swarm_params = zeros(npars,nswarms); % best parameters per swarm (used as the swarm's optimum)
best_particle = zeros(1,nswarms); % best particle number per swarm

% Continuous monitoring (per iteration) to check algorithm performance
global_opt_OFval_it = zeros(1,max_niter+1); % global optimum objective function value
global_opt_params_it = zeros(npars,max_niter+1); % global optimum parameters
local_opt_OFval_it = zeros(1,nparticles,max_niter+1); % local optimum objective function value
local_opt_params_it = zeros(npars,nparticles,max_niter+1); % local optimum parameters
particle_cur_params_it = zeros(npars,nparticles,max_niter+1); % particle current parameters
particle_energies_it = zeros(nparticles,max_niter); % particle energy

% For each swarm: select best particle and save
for swarm = 1:nswarms
    best_swarm_OFval(swarm) = min(particle_best_OFvals(swarm_groups(:,swarm)));
    best_particle(swarm) = find(particle_best_OFvals(swarm_groups(:,swarm)) == best_swarm_OFval(swarm),1);
    best_swarm_params(:,swarm) = particle_best_params(:,swarm_groups(best_particle(swarm),swarm));
    
%     load(['P_PSO',num2str(parallel_id),'_Particle_',num2str(swarm_groups(best_particle(swarm),swarm)),'.mat'],'P')
%     save(['P_BP_',datafile,'_swarm_',num2str(swarm),'.mat'],'P')
end

% Initialize particle velocities
velocities = zeros(npars,nparticles); % start with v = 0
norm_velocities = zeros(npars,nparticles); % normalized velocities for energy calculation
conv_arr = 1./(par_bounds(:,2)-par_bounds(:,1)); % conversion array for parameter normalization
conv_arr(isinf(conv_arr))=1; % change Inf to 1 for parameters which are not used in optimization to prevent error
particle_energies_par = ones(npars,nparticles); % particles energies per parameter

it = 0; % iteration 0
% best_swarm_OFval_it(it+1,:) = best_swarm_OFval;

global_opt_OFval_it(it+1) = min(particle_best_OFvals);
global_opt_params_it(:,it+1) = particle_best_params(:,find(particle_best_OFvals == min(particle_best_OFvals),1));
local_opt_OFval_it(1,:,it+1) = particle_best_OFvals; % local optimum objective function value
local_opt_params_it(:,:,it+1) = particle_best_params; % local optimum parameters
particle_cur_params_it(:,:,it+1) = particle_cur_params; % particle current parameters

while max(max(particle_energies_par)) >= 0.0001 && it < max_niter % Stop when particle energy is too low or maximum number of iterations is reached
    
    % Increase iteration counter
    it = it + 1;
    
    % Re-assign subswarms every 20 iterations (characteristic for 'dynamic' MS-PSO)
    if rem(it,20) == 1 && it ~= 1
        %disp('New random swarms formed')
        swarm_groups = reshape(randperm(nparticles),[nparticles/nswarms,nswarms]);
        
        % For each swarm: select best particle using the current positions
        for swarm = 1:nswarms
            best_swarm_OFval(swarm) = min(particle_cur_OFvals(swarm_groups(:,swarm)));
            best_particle(swarm) = find(particle_cur_OFvals(swarm_groups(:,swarm)) == best_swarm_OFval(swarm),1);
            best_swarm_params(:,swarm) = particle_cur_params(:,swarm_groups(best_particle(swarm),swarm));
        end
    end
    
    %disp(['Iteration: ',num2str(it)])
    %disp(['Best objective function value: ',num2str(min(particle_best_OFvals))])
    
    for swarm = 1:nswarms
        
        %disp(['Swarm ',num2str(swarm)])

        for particle = 1:nparticles/nswarms
            Z1 = rand(npars,1);
            Z2 = rand(npars,1);
            
            % Calculate velocities and normalized velocities
            norm_v = 0.729.*norm_velocities(:,swarm_groups(particle,swarm)) + 1.49445.*Z1.*(best_swarm_params(:,swarm)-particle_cur_params(:,swarm_groups(particle,swarm))).*conv_arr + 1.49445.*Z2.*(particle_best_params(:,swarm_groups(particle,swarm))-particle_cur_params(:,swarm_groups(particle,swarm))).*conv_arr;
            % Limit normalized velocity to 0.25 to prevent that too many particles fly outside of the boundaries
            norm_velocities(:,swarm_groups(particle,swarm)) = sign(norm_v).*min(abs(norm_v),0.25);
            
            % Calculate velocity
            velocities(:,swarm_groups(particle,swarm)) = norm_velocities(:,swarm_groups(particle,swarm))./conv_arr;
            
            % Calculate particle energies
            particle_energies_par(:,swarm_groups(particle,swarm)) = norm_velocities(:,swarm_groups(particle,swarm)).^2 + (conv_arr.*(particle_cur_params(:,swarm_groups(particle,swarm))-particle_best_params(:,swarm_groups(particle,swarm)))).^2 + (conv_arr.*(particle_cur_params(:,swarm_groups(particle,swarm))-best_swarm_params(:,swarm))).^2;
            
            %disp(['Particle ',num2str(particle),' energy: ',num2str(particle_energy)])

            % Update particle positions
            cur_pos = particle_cur_params(:,swarm_groups(particle,swarm));
            for step_np = 1:nsteps_np
                outb = 0; % out of boundary index, 0 = inside, 1 = outside
                particle_cur_params(:,swarm_groups(particle,swarm)) = cur_pos + step_np*(velocities(:,swarm_groups(particle,swarm))/nsteps_np);
                
                % Boundaries used for DMS-PSO (can be different than used during Monte Carlo)
                for ipar = 1:length(par_set)
                    if outb == 1
                        break
                    end
                    parname = par_set{ipar};
                
                    switch parname
                        % q0: [1,20] L/min
                        case 'q0'
                            if particle_cur_params(1,swarm_groups(particle,swarm)) < 1/(60*1000) || particle_cur_params(1,swarm_groups(particle,swarm)) > 20/(60*1000)
                                outb = 1;
                            end

                        % p0: [50,150] mmHg
                        case 'p0'
                            if particle_cur_params(2,swarm_groups(particle,swarm)) < 50*133.322 || particle_cur_params(2,swarm_groups(particle,swarm)) > 150*133.322
                                outb = 1;
                            end

                        % kSy: [1,1000] %
                        case 'k_Sy'
                            if particle_cur_params(3,swarm_groups(particle,swarm)) < 0.01 || particle_cur_params(3,swarm_groups(particle,swarm)) > 10
                                outb = 1;
                            end

                        % p0Sy: [1,1000] %
                        case 'p0_Sy'
                            if particle_cur_params(4,swarm_groups(particle,swarm)) < 0.01 || particle_cur_params(4,swarm_groups(particle,swarm)) > 10
                                outb = 1;
                            end

                        % A0Sy: [1,1000] %
                        case 'A0_Sy'
                            if particle_cur_params(5,swarm_groups(particle,swarm)) < 0.01 || particle_cur_params(5,swarm_groups(particle,swarm)) > 10
                                outb = 1;
                            end

                        % dTauAv: [-100,200] ms
                        case 'dTauAv'
                            if particle_cur_params(6,swarm_groups(particle,swarm)) < -0.100 || particle_cur_params(6,swarm_groups(particle,swarm)) > 0.200
                                outb = 1;
                            end

                        % dT Lv (18 patches)
                        case 'dT_Lv'
                        % (1) dT LVfw: [-60,200] ms
                            for seg = 1:12
                                if particle_cur_params(seg+6,swarm_groups(particle,swarm)) < -0.060 || particle_cur_params(seg+6,swarm_groups(particle,swarm)) > 0.200
                                    outb = 1;
                                end
                            end
                        % (2) dT septum: [-60,120] ms
                            for seg = 13:18
                                if particle_cur_params(seg+6,swarm_groups(particle,swarm)) < -0.060 || particle_cur_params(seg+6,swarm_groups(particle,swarm)) > 0.120
                                    outb = 1;
                                end
                            end
                            
                        case 'dT_g_Lv'
                            if particle_cur_params(7,swarm_groups(particle,swarm)) < -0.060 || particle_cur_params(7,swarm_groups(particle,swarm)) > 0.200
                                outb = 1;
                            end
        
                        % dT La: [-30,150] ms
                        case 'dT_La'
                            if particle_cur_params(25,swarm_groups(particle,swarm)) < -0.030 || particle_cur_params(25,swarm_groups(particle,swarm)) > 0.150
                                outb = 1;
                            end
                    
                        % VWall Lv (18 patches)
                        case 'VWall_Lv'
                        % (1) VWall LVfw: [20,300] %
                            for seg = 1:12
                                if particle_cur_params(seg+25,swarm_groups(particle,swarm)) < 0.2*8.0007e-6 || particle_cur_params(seg+25,swarm_groups(particle,swarm)) > 3*8.0007e-6
                                    outb = 1;
                                end
                            end
                        % (2) VWall septum: [20,300] %
                            for seg = 13:18
                                if particle_cur_params(seg+25,swarm_groups(particle,swarm)) < 0.2*5.3926e-6 || particle_cur_params(seg+25,swarm_groups(particle,swarm)) > 3*5.3926e-6
                                    outb = 1;
                                end
                            end
                        % VWall Lv (global): [20,300] %    
                        case 'VWall_g_Lv'
                            if particle_cur_params(26,swarm_groups(particle,swarm)) < 0.2 || particle_cur_params(26,swarm_groups(particle,swarm)) > 3
                                outb = 1;
                            end

                        % VWall La: [20,300] %
                        case 'VWall_La'
                            if particle_cur_params(44,swarm_groups(particle,swarm)) < 0.2*1.5583e-5 || particle_cur_params(44,swarm_groups(particle,swarm)) > 3*1.5583e-5
                                outb = 1;
                            end

                        % AmRef Lv (18 patches)
                        case 'AmRef_Lv'
                        % (1) AmRef LVfw: [20,300] %
                            for seg = 1:12
                                if particle_cur_params(seg+44,swarm_groups(particle,swarm)) < 0.2*8.1711e-4 || particle_cur_params(seg+44,swarm_groups(particle,swarm)) > 3*8.1711e-4
                                    outb = 1;
                                end
                            end
                        % (2) AmRef septum: [20,300] %
                            for seg = 13:18
                                if particle_cur_params(seg+44,swarm_groups(particle,swarm)) < 0.2*8.1404e-4 || particle_cur_params(seg+44,swarm_groups(particle,swarm)) > 3*8.1404e-4
                                    outb = 1;
                                end
                            end
                        % AmRef Lv (global): [20,300] %    
                        case 'AmRef_g_Lv'
                            if particle_cur_params(45,swarm_groups(particle,swarm)) < 0.2 || particle_cur_params(45,swarm_groups(particle,swarm)) > 3
                                outb = 1;
                            end

                        % AmRef Rv: [20,300] %
                        case 'AmRef_Rv'
                            if particle_cur_params(63,swarm_groups(particle,swarm)) < 0.2*0.0130 || particle_cur_params(63,swarm_groups(particle,swarm)) > 3*0.0130
                                outb = 1;
                            end

                        % AmRef La: [20,300] %
                        case 'AmRef_La'
                            if particle_cur_params(64,swarm_groups(particle,swarm)) < 0.2*0.0069 || particle_cur_params(64,swarm_groups(particle,swarm)) > 3*0.0069
                                outb = 1;
                            end

                        % SfAct Lv (18 patches): [0,1000] kPa
                        case 'SfAct_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+64,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+64,swarm_groups(particle,swarm)) > 1e6
                                    outb = 1;
                                end
                            end
                        % SfAct Lv (global): [0,1000] kPa    
                        case 'SfAct_g_Lv'
                            if particle_cur_params(65,swarm_groups(particle,swarm)) < 0 || particle_cur_params(65,swarm_groups(particle,swarm)) > 1e6
                                outb = 1;
                            end

                        % SfAct Rv: [0,1000] kPa
                        case 'SfAct_Rv'
                            if particle_cur_params(83,swarm_groups(particle,swarm)) < 0 || particle_cur_params(83,swarm_groups(particle,swarm)) > 1e6
                                outb = 1;
                            end

                        % SfAct La: [0,1000] kPa
                        case 'SfAct_La'
                            if particle_cur_params(84,swarm_groups(particle,swarm)) < 0 || particle_cur_params(84,swarm_groups(particle,swarm)) > 1e6
                                outb = 1;
                            end

                        % vMax Lv (18 patches): [0.5,50] microm/s
                        case 'vMax_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+84,swarm_groups(particle,swarm)) < 0.5 || particle_cur_params(seg+84,swarm_groups(particle,swarm)) > 50
                                    outb = 1;
                                end
                            end
                        % vMax Lv (global): [0.5,50] microm/s    
                        case 'vMax_g_Lv'
                            if particle_cur_params(85,swarm_groups(particle,swarm)) < 0.5 || particle_cur_params(85,swarm_groups(particle,swarm)) > 50
                                outb = 1;
                            end

                        % TR Lv (18 patches): [0,1]
                        case 'TR_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+102,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+102,swarm_groups(particle,swarm)) > 1
                                    outb = 1;
                                end
                            end
                        % TR Lv (global): [0,1]    
                        case 'TR_g_Lv'
                            if particle_cur_params(103,swarm_groups(particle,swarm)) < 0 || particle_cur_params(103,swarm_groups(particle,swarm)) > 1
                                outb = 1;
                            end

                        % TR Rv: [0,1]
                        case 'TR_Rv'
                            if particle_cur_params(121,swarm_groups(particle,swarm)) < 0 || particle_cur_params(121,swarm_groups(particle,swarm)) > 1
                                outb = 1;
                            end

                        % TR La: [0,1]
                        case 'TR_La'
                            if particle_cur_params(122,swarm_groups(particle,swarm)) < 0 || particle_cur_params(122,swarm_groups(particle,swarm)) > 1
                                outb = 1;
                            end

                        % TD Lv (18 patches): [0,1]
                        case 'TD_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+122,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+122,swarm_groups(particle,swarm)) > 1
                                    outb = 1;
                                end
                            end
                        % TD Lv (global): [0,1]    
                        case 'TD_g_Lv'
                            if particle_cur_params(123,swarm_groups(particle,swarm)) < 0 || particle_cur_params(123,swarm_groups(particle,swarm)) > 1
                                outb = 1;
                            end

                        % Ls0Pas Lv (18 patches): [0,3] microm
                        case 'Ls0Pas_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+140,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+140,swarm_groups(particle,swarm)) > 3
                                    outb = 1;
                                end
                            end
                        % Ls0Pas Lv (global): [0,3] microm    
                        case 'Ls0Pas_g_Lv'
                            if particle_cur_params(141,swarm_groups(particle,swarm)) < 0 || particle_cur_params(141,swarm_groups(particle,swarm)) > 3
                                outb = 1;
                            end

                        % Ls0Pas La: [0,3] microm
                        case 'Ls0Pas_La'
                            if particle_cur_params(159,swarm_groups(particle,swarm)) < 0 || particle_cur_params(159,swarm_groups(particle,swarm)) > 3
                                outb = 1;
                            end

                        % SfPas Lv (18 patches): [0,1000] kPa
                        case 'SfPas_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+159,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+159,swarm_groups(particle,swarm)) > 1e6
                                    outb = 1;
                                end
                            end
                        % SfPas Lv (global): [0,1000] kPa    
                        case 'SfPas_g_Lv'
                            if particle_cur_params(160,swarm_groups(particle,swarm)) < 0 || particle_cur_params(160,swarm_groups(particle,swarm)) > 1e6/2.2311e4
                                outb = 1;
                            end

                        % k1 Lv (18 patches): [0,100]
                        case 'k1_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+177,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+177,swarm_groups(particle,swarm)) > 100
                                    outb = 1;
                                end
                            end
                        % k1 Lv (global): [0,100]    
                        case 'k1_g_Lv'
                            if particle_cur_params(178,swarm_groups(particle,swarm)) < 0 || particle_cur_params(178,swarm_groups(particle,swarm)) > 100
                                outb = 1;
                            end

                        % k1 La: [0,100]
                        case 'k1_La'
                            if particle_cur_params(196,swarm_groups(particle,swarm)) < 0 || particle_cur_params(196,swarm_groups(particle,swarm)) > 100
                                outb = 1;
                            end

                        % dLsPas Lv (18 patches): [0.01,10]
                        case 'dLsPas_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+196,swarm_groups(particle,swarm)) < 0.01 || particle_cur_params(seg+196,swarm_groups(particle,swarm)) > 10
                                    outb = 1;
                                end
                            end
                        % dLsPas Lv (global): [0.01,10]    
                        case 'dLsPas_g_Lv'
                            if particle_cur_params(197,swarm_groups(particle,swarm)) < 0.01 || particle_cur_params(197,swarm_groups(particle,swarm)) > 10
                                outb = 1;
                            end

                        % LenSeriesElement LV global: [0,0.2]
                        case 'LenSeriesElement_g_Lv'
                            if particle_cur_params(215,swarm_groups(particle,swarm)) < 0 || particle_cur_params(215,swarm_groups(particle,swarm)) > 0.2
                                outb = 1;
                            end

                        % ADO Lv (18 patches): [0,cycletime*2] ms
                        case 'ADO_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+215,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+215,swarm_groups(particle,swarm)) > measurement.ctmean*2
                                    outb = 1;
                                end
                            end
                        % ADO Lv (global): [0,cycletime*2] ms    
                        case 'ADO_g_Lv'
                            if particle_cur_params(216,swarm_groups(particle,swarm)) < 0 || particle_cur_params(216,swarm_groups(particle,swarm)) > measurement.ctmean*2
                                outb = 1;
                            end

                        % ADO La: [0,cycletime*2] ms
                        case 'ADO_La'
                            if particle_cur_params(234,swarm_groups(particle,swarm)) < 0 || particle_cur_params(234,swarm_groups(particle,swarm)) > measurement.ctmean*2
                                outb = 1;
                            end

                        % LDAD Lv (18 patches): [0,2]
                        case 'LDAD_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+234,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+234,swarm_groups(particle,swarm)) > 2
                                    outb = 1;
                                end
                            end
                        % LDAD Lv (global): [0,2]    
                        case 'LDAD_g_Lv'
                            if particle_cur_params(235,swarm_groups(particle,swarm)) < 0 || particle_cur_params(235,swarm_groups(particle,swarm)) > 2
                                outb = 1;
                            end

                        % LDCI Lv (18 patches): [0,50]
                        case 'LDCI_Lv'
                            for seg = 1:18
                                if particle_cur_params(seg+252,swarm_groups(particle,swarm)) < 0 || particle_cur_params(seg+252,swarm_groups(particle,swarm)) > 50
                                    outb = 1;
                                end
                            end
                        % LDCI Lv (global): [0,50]    
                        case 'LDCI_g_Lv'
                            if particle_cur_params(253,swarm_groups(particle,swarm)) < 0 || particle_cur_params(253,swarm_groups(particle,swarm)) > 50
                                outb = 1;
                            end
                    end
                end
 
                if outb ~= 1
                    % Set parameter values
                    P = setSim(['P_PSO',num2str(parallel_id),'_Particle',num2str(swarm_groups(particle,swarm)),'.mat'],par_set,particle_cur_params(:,swarm_groups(particle,swarm)));
                    
                    try
                        RunCPP1
                    catch
                        %disp('Particle did not run successfully')
                    end
                else
                end
            end
            if outb ~= 1
                % Calculate objective function value
                try
                    if any(isnan(P.Valve.q(:)))
                        error('NaN found in q')
                    end
                    
                    if min(P.Cavity.V(:,7).*10^6) <= 0 || min(P.Cavity.V(:,8).*10^6) <= 0
                        error('Volume criteria not met')
                    end
                    
                    if min(P.Node.p(:,8)./133.322) <= 0 || min(P.Node.p(:,7)./133.322) <= 0
                        error('Pressure criteria not met')
                    end

                    [OFval_i,~] = CalculateObjectiveFunction_18(P,measurement);
                    particle_cur_OFvals(swarm_groups(particle,swarm)) = OFval_i;
                    
                    save(['P_PSO',num2str(parallel_id),'_Particle',num2str(swarm_groups(particle,swarm)),'.mat'],'P')

                    if OFval_i < particle_best_OFvals(swarm_groups(particle,swarm))
                        particle_best_OFvals(swarm_groups(particle,swarm)) = OFval_i;
                        particle_best_params(:,swarm_groups(particle,swarm)) = particle_cur_params(:,swarm_groups(particle,swarm));
                        if OFval_i == min(particle_best_OFvals)
                            save(['P_BestParticle_',datafile,'_',num2str(rep_id),'.mat'],'P')
                        end
                    end

                catch
                    particle_cur_OFvals(swarm_groups(particle,swarm)) = Inf;
                    CA = CircAdapt();
                    CA.buildTriSeg();
                    CA=setP(CA,P);
                end
            else
                particle_cur_OFvals(swarm_groups(particle,swarm)) = Inf;
            end
        end
    
        if min(particle_cur_OFvals(swarm_groups(:,swarm))) < best_swarm_OFval(swarm)
            % update swarm best position
            best_particle(swarm) = find(particle_cur_OFvals(swarm_groups(:,swarm)) == min(particle_cur_OFvals(swarm_groups(:,swarm))),1);
            best_swarm_params(:,swarm) = particle_cur_params(:,swarm_groups(best_particle(swarm),swarm));
            best_swarm_OFval(swarm) = min(particle_cur_OFvals(swarm_groups(:,swarm)));
        end
    end
    global_opt_OFval_it(it+1) = min(particle_best_OFvals);
    global_opt_params_it(:,it+1) = particle_best_params(:,find(particle_best_OFvals == min(particle_best_OFvals),1));
    local_opt_OFval_it(1,:,it+1) = particle_best_OFvals; % local optimum objective function value
    local_opt_params_it(:,:,it+1) = particle_best_params; % local optimum parameters
    particle_cur_params_it(:,:,it+1) = particle_cur_params; % particle current parameters
    particle_energies_it(:,it) = sum(particle_energies_par,1)';
end

global_opt_params = global_opt_params_it(:,it+1);
global_opt_OFval = global_opt_OFval_it(it+1);

% Save output of interest to file
save(['Results_PSO_',datafile,'_',num2str(rep_id),'.mat'],'global_opt_OFval_it','global_opt_params_it','local_opt_OFval_it','local_opt_params_it','particle_cur_params_it','particle_energies_it')
            