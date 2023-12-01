% mediation analysis -----------------------------------------------------------
% peer environment of ABCD -----------------------------------------------------
% RSFCs - Behaviors ------------------------------------------------------------
clear all; clc;
cd H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/mediations

% load data --------------------------------------------------------------------
mediations_rsfcs = readtable('mediation_rsfc.csv');
vars = mediations_rsfcs.Properties.VariableNames(2:120)';
mediations = table2array(mediations_rsfcs(:,2:120));

% PFI - cortical RSFCs
PFI_ad_ca = PFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_ad_ngd_ca')));
PFI_ad_sa = PFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_ad_ngd_sa')));
PFI_cgc_dla = PFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_cgc_ngd_dla')));
PFI_ca_smh = PFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_ca_ngd_smh')));
PFI_ca_smm = PFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_ca_ngd_smm')));
PFI_ca_vs = PFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_ca_ngd_vs')));
PFI_dt_fo = PFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_dt_ngd_fo')));
PFI_dt_rspltp = PFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_dt_ngd_rspltp')));

% DFI - cortical RSFCs
DFI_dt_dt = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_dt_ngd_dt')));
DFI_dla_dla = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_dla_ngd_dla')));
DFI_smh_smh = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_smh_ngd_smh')));
DFI_cgc_dt = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_cgc_ngd_dt')));
DFI_cgc_vs = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_cgc_ngd_vs')));
DFI_dt_dla = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_dt_ngd_dla')));
DFI_dla_vta = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_dla_ngd_vta')));
DFI_smh_smm = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_c_ngd_smh_ngd_smm')));

% DFI - cortico-subcortical RSFCs
DFI_au_cde = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_au_scs_cde')));
DFI_au_hp = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_au_scs_hp')));
DFI_au_aa = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_au_scs_aa')));
DFI_cerc_thp = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_cerc_scs_thp')));
DFI_cerc_pt = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_cerc_scs_pt')));
DFI_cerc_hp = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_cerc_scs_hp')));
DFI_cerc_ag = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_cerc_scs_ag')));
DFI_copa_crcx = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_copa_scs_crcx')));
DFI_copa_cde = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_copa_scs_cde')));
DFI_copa_vtdc = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_copa_scs_vtdc')));
DFI_df_pl = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_df_scs_pl')));
DFI_df_ag = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_df_scs_ag')));
DFI_df_aa = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_df_scs_aa')));
DFI_dsa_hp = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_dsa_scs_hp')));
DFI_dsa_ag = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_dsa_scs_ag')));
DFI_dsa_vtdc = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_dsa_scs_vtdc')));
DFI_fopa_crcx = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_fopa_scs_crcx')));
DFI_fopa_thp = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_fopa_scs_thp')));
DFI_fopa_ag = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_fopa_scs_ag')));
DFI_rst_aa = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_rst_scs_aa')));
DFI_smh_crcx = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smh_scs_crcx')));
DFI_smh_cde = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smh_scs_cde')));
DFI_smh_pt = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smh_scs_pt')));
DFI_smh_pl = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smh_scs_pl')));
DFI_smh_aa = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smh_scs_aa')));
DFI_smm_pl = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smm_scs_pl')));
DFI_smm_hp = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smm_scs_hp')));
DFI_smm_ag = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smm_scs_ag')));
DFI_sa_hp = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_sa_scs_hp')));
DFI_sa_aa = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_sa_scs_aa')));
DFI_sa_vtdc = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_sa_scs_vtdc')));
DFI_vta_cde = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_vta_scs_cde')));
DFI_cerc_bs = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_cerc_scs_bs')));
DFI_dsa_bs = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_dsa_scs_bs')));
DFI_rst_bs = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_rst_scs_bs')));
DFI_smh_bs = DFI_rsfcs(mediations, find(strcmp(vars, 'rsfmri_cor_ngd_smh_scs_bs')));

% output ------------------------------------------------------------------
% save results
save meadiations_rsfcs_results.mat

% PFI
csvwrite('PFI_ad_ca.csv', PFI_ad_ca)
csvwrite('PFI_ad_sa.csv', PFI_ad_sa)
csvwrite('PFI_cgc_dla.csv', PFI_cgc_dla)
csvwrite('PFI_ca_smh.csv', PFI_ca_smh)
csvwrite('PFI_ca_smm.csv', PFI_ca_smm)
csvwrite('PFI_ca_vs.csv', PFI_ca_vs)
csvwrite('PFI_dt_fo.csv', PFI_dt_fo)
csvwrite('PFI_dt_rspltp.csv', PFI_dt_rspltp)

% DFI
csvwrite('DFI_dt_dt.csv', DFI_dt_dt)
csvwrite('DFI_dla_dla.csv', DFI_dla_dla)
csvwrite('DFI_smh_smh.csv', DFI_smh_smh)
csvwrite('DFI_cgc_dt.csv', DFI_cgc_dt)
csvwrite('DFI_cgc_vs.csv', DFI_cgc_vs)
csvwrite('DFI_dt_dla.csv', DFI_dt_dla)
csvwrite('DFI_dla_vta.csv', DFI_dla_vta)
csvwrite('DFI_smh_smm.csv', DFI_smh_smm)
% DFI cortico-subcortical FCs
csvwrite('DFI_au_cde.csv', DFI_au_cde)
csvwrite('DFI_au_hp.csv', DFI_au_hp)
csvwrite('DFI_au_aa.csv', DFI_au_aa)
csvwrite('DFI_cerc_thp.csv', DFI_cerc_thp)
csvwrite('DFI_cerc_pt.csv', DFI_cerc_pt)
csvwrite('DFI_cerc_hp.csv', DFI_cerc_hp)
csvwrite('DFI_cerc_ag.csv', DFI_cerc_ag)
csvwrite('DFI_copa_crcx.csv', DFI_copa_crcx)
csvwrite('DFI_copa_cde.csv', DFI_copa_cde)
csvwrite('DFI_copa_vtdc.csv', DFI_copa_vtdc)
csvwrite('DFI_df_pl.csv', DFI_df_pl)
csvwrite('DFI_df_ag.csv', DFI_df_ag)
csvwrite('DFI_df_aa.csv', DFI_df_aa)
csvwrite('DFI_dsa_hp.csv', DFI_dsa_hp)
csvwrite('DFI_dsa_ag.csv', DFI_dsa_ag)
csvwrite('DFI_dsa_vtdc.csv', DFI_dsa_vtdc)
csvwrite('DFI_fopa_crcx.csv', DFI_fopa_crcx)
csvwrite('DFI_fopa_thp.csv', DFI_fopa_thp)
csvwrite('DFI_fopa_ag.csv', DFI_fopa_ag)
csvwrite('DFI_rst_aa.csv', DFI_rst_aa)
csvwrite('DFI_smh_crcx.csv', DFI_smh_crcx)
csvwrite('DFI_smh_cde.csv', DFI_smh_cde)
csvwrite('DFI_smh_pt.csv', DFI_smh_pt)
csvwrite('DFI_smh_pl.csv', DFI_smh_pl)
csvwrite('DFI_smh_aa.csv', DFI_smh_aa)
csvwrite('DFI_smm_pl.csv', DFI_smm_pl)
csvwrite('DFI_smm_hp.csv', DFI_smm_hp)
csvwrite('DFI_smm_ag.csv', DFI_smm_ag)
csvwrite('DFI_sa_hp.csv', DFI_sa_hp)
csvwrite('DFI_sa_aa.csv', DFI_sa_aa)
csvwrite('DFI_sa_vtdc.csv', DFI_sa_vtdc)
csvwrite('DFI_vta_cde.csv', DFI_vta_cde)
csvwrite('DFI_cerc_bs.csv', DFI_cerc_bs)
csvwrite('DFI_dsa_bs.csv', DFI_dsa_bs)
csvwrite('DFI_rst_bs.csv', DFI_rst_bs)
csvwrite('DFI_smh_bs.csv', DFI_smh_bs)