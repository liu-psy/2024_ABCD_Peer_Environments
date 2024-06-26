function [results] = PFI_rsfcs_brain(data, Y_index)
    % covariates
    cov = data(:,1:8);
    
    PFI = zscore(data(:,9));
    brain = zscore(data(:,Y_index));

    R1 = []; R2 = []; R3 = [];
    Path1 = []; Path2 = []; Path3 = []; 
    Stats = {};
    
    all_behavior = data(:, 11:67);
  
    for i = 1:57
        behavior = zscore(all_behavior(:,i));
        [paths1, stats1] = mediation(PFI, behavior, brain, 'boot', 'bootsamples', 10000, 'covs', cov);
        R1(i,:) = stats1.p;
        Path1(i,:) = paths1;
        CIs(i,:) = [stats1.ci(:,:,1) stats1.ci(:,:,2)];
        
        i
    end
results = [R1 Path1 CIs]; 