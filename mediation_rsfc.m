% mediation analysis -----------------------------------------------------------
% peer environment of ABCD -----------------------------------------------------
% RSFCs - Behaviors ------------------------------------------------------------
clear all; clc;
cd H:/ABCD/Release4.0/Package_1194636/results/peer_environments/mediations

% load data --------------------------------------------------------------------
mediations_rsfcs = readtable('mediation_rsfc.csv');
vars = mediations_rsfcs.Properties.VariableNames(2:width(mediations_rsfcs))';
mediations = table2array(mediations_rsfcs(:,2:width(mediations_rsfcs)));

PFI_RSFCs_fdr = readtable('PFI_RSFC_fdr.txt');                             % significant RSFCs for PFI
DFI_RSFCs_fdr = readtable('DFI_RSFC_fdr.txt');                             % significant RSFCs for DFI

% PFI - 8+4 RSFCs
PFI_rsfc_results = {};
for i = 1:length(PFI_RSFCs_fdr.x)
    i
    mediation_results = PFI_rsfcs(mediations, find(strcmp(vars, PFI_RSFCs_fdr.x{i})));
    PFI_rsfc_results{i} = mediation_results;
    % get FC
    fc = strsplit(PFI_RSFCs_fdr.x{i}, '_');
    file_names = strcat('PFI_', fc{4}, '_', fc{6}, '.csv');
    
    % output
    csvwrite(file_names, mediation_results)
end


% DFI - 8+36 RSFCs
DFI_rsfc_results = {};
for i = 1:length(DFI_RSFCs_fdr.x)
    i
    mediation_results = DFI_rsfcs(mediations, find(strcmp(vars, DFI_RSFCs_fdr.x{i})));
    DFI_rsfc_results{i} = mediation_results;
    % get FC
    fc = strsplit(DFI_RSFCs_fdr.x{i}, '_');
    file_names = strcat('DFI_', fc{4}, '_', fc{6}, '.csv');
    
    % output
    csvwrite(file_names, mediation_results)
end

save meadiations_rsfcs_results.mat