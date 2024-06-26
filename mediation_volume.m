% mediation analysis -----------------------------------------------------------
% peer environment of ABCD -----------------------------------------------------
% Brain strcutre - Behaviors ---------------------------------------------------
clear all; clc;
cd H:/ABCD/Release4.0/Package_1194636/results/peer_environments/mediations

% load data
mediations_volumes = readtable('mediation_volume.csv');
vars = mediations_volumes.Properties.VariableNames(2:width(mediations_volumes))';
mediations = table2array(mediations_volumes(:,2:width(mediations_volumes)));

PFI_volume_fdr = readtable('PFI_volume_fdr.txt');                          % significant volumes for PFI
DFI_volume_fdr = readtable('DFI_volume_fdr.txt');                          % significant volumes for DFI

% mediation --------------------------------------------------------------------
% PFI (10 regions)
PFI_volume_results = {};
for i = 1:length(PFI_volume_fdr.x) 
    i
    mediation_results = PFI_volumes(mediations, find(strcmp(vars, PFI_volume_fdr.x{i})));
    PFI_volume_results{i} = mediation_results;
    % get FC
    volume = strsplit(PFI_volume_fdr.x{i}, '_');
    file_names = strcat('PFI_', volume{4}, '.csv');
    
    % output
    csvwrite(file_names, mediation_results)
end

% DFI (one regions)
DFI_volume_results = {};
for i = 1:length(DFI_volume_fdr.x) 
    i
    mediation_results = PFI_volumes(mediations, find(strcmp(vars, DFI_volume_fdr.x{i})));
    PFI_volume_results{i} = mediation_results;
    % get FC
    volume = strsplit(DFI_volume_fdr.x{i}, '_');
    file_names = strcat('DFI_', volume{4}, '.csv');
    
    % output
    csvwrite(file_names, mediation_results)
end

% save results
save mediation_volume_results.mat