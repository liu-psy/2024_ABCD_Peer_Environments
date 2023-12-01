% This script is used to plot brain for Figure2

clear; clc
cd H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments

% load data (68*2)
dk_table = readmatrix('dk_table.txt');

% Map parcellated data to the surface (volume)
PFI = parcel_to_surface(dk_table(:, 1), 'aparc_fsa5');
DFI = parcel_to_surface(dk_table(:, 2), 'aparc_fsa5');

% Project the results on the surface brain
p1 = figure,
    plot_cortical(PFI, 'surface_name', 'fsa5', 'color_range', ...
    [-4.5 4.5], 'cmap', 'RdBu_r');
exportgraphics(p1, 'figures/cortical_PFI.tif', 'Resolution', 600);

p2 = figure,
    plot_cortical(DFI, 'surface_name', 'fsa5', 'color_range', ...
    [-4.5 4.5], 'cmap', 'RdBu_r');
exportgraphics(p2, 'figures/cortical_DFI.tif', 'Resolution', 600);