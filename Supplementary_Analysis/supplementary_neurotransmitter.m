% Spatial correlations between t-maps of brain structure-peer environments associations and neurotransmitter density
% Supplemetary analysis
% Toolbox: JuSpace

clear all; clc;

cd 'H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/supplementary_analysis'
% create t maps
ref_num = [1001;1002;1003;1004;1006;1007;1008;1009;1010;1011;1012;1013;1014;1015;1016;1017;1018;...
    1019;1020;1021;1022;1023;1024;1025;1026;1027;1028;1029;1030;1031;1032;1033;1034;1035;...
    2001;2002;2003;2004;2006;2007;2008;2009;2010;2011;2012;2013;2014;2015;2016;2017;2018;...
    2019;2020;2021;2022;2023;2024;2025;2026;2027;2028;2029;2030;2031;2032;2033;2034;2035];

APARC_atlas = 'aparc.nii';

% load t maps from LMM models
tmap_vol_DFI = importdata('tmap_vol_DFI.txt');
tmap_vol_PFI = importdata('tmap_vol_PFI.txt');
tmap_thick_DFI = importdata('tmap_thick_DFI.txt');
tmap_thick_PFI = importdata('tmap_thick_PFI.txt');
tmap_area_DFI = importdata('tmap_area_DFI.txt');
tmap_area_PFI = importdata('tmap_area_PFI.txt');

% ROT to Nifti
ROI2nifti(tmap_vol_DFI, APARC_atlas, ref_num, 'vol_DFI');
ROI2nifti(tmap_vol_PFI, APARC_atlas, ref_num, 'vol_PFI');
ROI2nifti(tmap_thick_DFI, APARC_atlas, ref_num, 'thick_DFI');
ROI2nifti(tmap_thick_PFI, APARC_atlas, ref_num, 'thick_PFI');
ROI2nifti(tmap_area_DFI, APARC_atlas, ref_num, 'area_DFI');
ROI2nifti(tmap_area_PFI, APARC_atlas, ref_num, 'area_PFI');

% command for JuSpace -----------------------------------------------------
group1 = {
    'area_DFI.nii'
    'thick_DFI.nii'
    'vol_DFI.nii'
    'area_PFI.nii'
    'thick_PFI.nii'
    'vol_PFI.nii'
    };
group2 = {};
PET = {
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\5HT1a_WAY_HC36.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\5HT1b_P943_HC22.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\5HT2a_ALT_HC19.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\SERT_DASB_HC30.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\D1_SCH23390_c11.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\D2_RACLOPRIDE_c11.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\DAT_DATSPECT.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\FDOPA_f18.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\GABAa_FLUMAZENIL_c11.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\NAT_MRB_c11.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\VAChT_feobv_hc18_aghourian.nii'
        'E:\Matlab\Packages\JuSpace\JuSpace_v1.4\PETatlas\mGluR5_abp_hc73_smart.nii'
        };
% set atlas (dk)
ATLAS = 'E:/Matlab/Packages/JuSpace/JuSpace_v1.4/atlas/aparc.nii';
% each file in list 1
OPTIONS = [4; 1; 0; 0; 1];
% permutaion
Nperm = 1000;

% compute
[res,p_all,stats,data,D1,D2,data_PET,Resh,T1] = compute_DomainGauges(group1, group2, PET, ATLAS, OPTIONS, 'PET_result');
% compute p_exact
[p_exact,dist_rand] = compute_exact_spatial_pvalue(D1,data_PET, ATLAS, res, Nperm, OPTIONS, PET, T1);
% results (add p_exact)
Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];

% output results
writecell(Resh, 'Resh.csv')
% save results
save Supplemetary_Results_Neurotransmitter.mat ...
    Resh Nperm stats D1 data data_PET OPTIONS dist_rand