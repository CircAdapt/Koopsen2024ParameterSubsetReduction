% ---CalculateMaxDiaphony.m---
%
% Calculate maximum diaphony for all parameter groups
%
% Koopsen T.
% Last modified: 11/28/2023

% Patient datasets
datanames = {'MARC001','MARC008','MARC009','MARC022','MARC063','MARC072','MARC076','DEFIMI33','DEFIMI40','DEFIMI41','DEFIMI44','DEFIMI57','DEFIMI84'};

input = sobolset(103,'Skip',1000); % input Sobol simulations (has to match definition in SobolSamples.m (line 179))

% Define parameter groups based on the current parameter subset (to be modified when using a different parameter subset)
pargroups = {1,2,3,4,5,6,7,8:13,14,15:20,21,22:27,28,29,30:35,36,37,38:43,44:49,50,51,52:57,58:63,64,65:70,71:76,77,78:83,84,85:90,91,92:97,98:103};

% Determine the total number of parameters in the current subset
npars = 0;
for i = 1:length(pargroups)
    npars = npars + length(pargroups{i});
end

Nbest = 2*npars; % number of best simulations to be determined (twice the total number of parameters)
ncomp = 15; % number of objective function value components, diaphony is calculated for each component and for the sum of all 15 components

for d = 1:length(datanames)
    load(['OF_',datanames{d},'.mat'],'all_OF')
    eval(['all_OF',num2str(d),'=all_OF;'])
end

% Bootstrapping can be applied to determine confidence intervals, but not
% used here. Therefore, nboot is set to 1.
nboot = 1;

diaphony_per_patient = zeros(length(pargroups),length(datanames),ncomp+1,nboot);
diaphony_overall = zeros(length(pargroups),ncomp+1,nboot);
diaphony_comps = zeros(length(pargroups),ncomp+1);

for iboot = 1:nboot
    
    % In case of bootstrapping: use this line
    %samples = randi(3000000,3000000,1);
    
    % No bootstrapping
    samples = 1:3000000;

    for comp = 1:ncomp+1
        for d = 1:length(datanames)
            if comp<=ncomp
                % Individual component of objective function (15)
                eval(['all_OF_OFcomp=all_OF',num2str(d),'(:,:,comp);'])
            else
                % Sum of all components of objective function (1)
                eval(['all_OF_OFcomp=sum(all_OF',num2str(d),',3);'])
            end
            all_OF_samples = all_OF_OFcomp(samples);

            max_diaphony_per_pargroup = zeros(length(pargroups),1);
    
            [~,I] = mink(all_OF_samples,Nbest);
            index = samples(I);
    
            % Calculate parameter diaphony
            diaphony = CalculateDiaphony(index,input);
            
            % Calculate maximum diaphony per parameter group
            for n = 1:length(pargroups)
                max_diaphony_per_pargroup(n) = max(diaphony(pargroups{n}));
            end
            diaphony_per_patient(:,d,comp,iboot) = max_diaphony_per_pargroup;

            clear all_OF_samples
        end

        diaphony_overall(:,comp,iboot) = max(diaphony_per_patient(:,:,comp,iboot),[],2);
        mean_diaphony_overall = mean(diaphony_overall,3);
        diaphony_comps(:,comp) = mean_diaphony_overall(:,comp);
    end  
    %ci_diaphony_overall = prctile(diaphony_overall,[5 50 95],2);
    %std_diaphony_overall = std(diaphony_overall,0,2);
end
