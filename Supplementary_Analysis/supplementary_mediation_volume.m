% mediation analysis -----------------------------------------------------------
% peer environment of ABCD -----------------------------------------------------
% Brain strcutre - Behaviors ---------------------------------------------------
clear all; clc;
cd H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/supplementary_analysis/mediations

% load data
mediations_volumes = readtable('mediation_volume.csv');
vars = mediations_volumes.Properties.VariableNames(2:82)';
mediations = table2array(mediations_volumes(:,2:82));

% mediation --------------------------------------------------------------------
% PFI (4 regions)
PFI_ifpllh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_ifpllh')));
PFI_precnlh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_precnlh')));

PFI_putamenlh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_scs_putamenlh')));
PFI_pallidumrh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_scs_pallidumrh')));

% DFI (one region)
DFI_locclh = DFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_locclh')));

% output files ------------------------------------------------------------
save mediation_volume_results.mat

% output files
csvwrite('PFI_ifpllh.csv', PFI_ifpllh)
csvwrite('PFI_precnlh.csv', PFI_precnlh)
csvwrite('PFI_putamenlh.csv', PFI_putamenlh)
csvwrite('PFI_pallidumrh.csv', PFI_pallidumrh)

csvwrite('DFI_locclh.csv', DFI_locclh)