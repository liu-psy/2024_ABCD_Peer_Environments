function writeNiftiRef(Data, AtlasPath, RefNum, ResultName)
V = spm_vol(AtlasPath); [Y,~] = spm_read_vols(V);
vlat = V;vlat.fname = [ResultName,'.nii'];
voltmplat = zeros(size(Y));
D = length(Data);
for i = 1:D
    voltmplat(find(Y == RefNum(i))) = Data(i);
end
vlat.dt = [16,0];
vlat.descrip = 'FSL5.0';
spm_write_vol(vlat,voltmplat);