% mediation analysis -----------------------------------------------------------
% peer environment of ABCD -----------------------------------------------------
% Brain strcutre - Behaviors ---------------------------------------------------
clear all; clc;
cd H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/mediations

% load data
mediations_volumes = readtable('mediation_volume.csv');
vars = mediations_volumes.Properties.VariableNames(2:80)';
mediations = table2array(mediations_volumes(:,2:80));

% mediation --------------------------------------------------------------------
% PFI (10 regions)
PFI_ifpllh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_ifpllh')));
PFI_paracnlh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_paracnlh')));
PFI_precnlh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_precnlh')));
PFI_sufrlh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_sufrlh')));
PFI_postcnrh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_postcnrh')));
PFI_insularh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_insularh')));

PFI_putamenlh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_scs_putamenlh')));
PFI_putamenrh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_scs_putamenrh')));
PFI_pallidumrh = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_scs_pallidumrh')));
PFI_aar = PFI_volumes(mediations, find(strcmp(vars, 'smri_vol_scs_aar')));

% DFI (one regions)
DFI_locclh = DFI_volumes(mediations, find(strcmp(vars, 'smri_vol_cdk_locclh')));

% output files ------------------------------------------------------------
save mediation_volume_results.mat

% output files
csvwrite('PFI_ifpllh.csv', PFI_ifpllh)
csvwrite('PFI_paracnlh.csv', PFI_paracnlh)
csvwrite('PFI_precnlh.csv', PFI_precnlh)
csvwrite('PFI_sufrlh.csv', PFI_sufrlh)
csvwrite('PFI_postcnrh.csv', PFI_postcnrh)
csvwrite('PFI_insularh.csv', PFI_insularh)
csvwrite('PFI_putamenlh.csv', PFI_putamenlh)
csvwrite('PFI_putamenrh.csv', PFI_putamenrh)
csvwrite('PFI_pallidumrh.csv', PFI_pallidumrh)
csvwrite('PFI_aar.csv', PFI_aar)

csvwrite('DFI_locclh.csv', DFI_locclh)
