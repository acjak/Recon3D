close all; clear;

addpath('npy_matlab_master/');
% Read reconstructed volume. Format: X, Y, Z, param. Parameters: gamma, mu,
% completeness
V = readNPY('/u/data/alcer/DFXRM_rec/Rec_test_2/grain_ang.npy');

% Plot shape distribution
V_shape = zeros(size(V,1), size(V,2), size(V,3));
V_th = zeros(size(V,1), size(V,2), size(V,3));
for ii =1:size(V,1)
    for jj = 1:size(V,2)
        for kk = 1:size(V,3)
            V_shape(ii,jj,kk) = V(ii,jj,kk,3);
            if V(ii,jj,kk,3) > 0.75
                V_th(ii,jj,kk) = V(ii,jj,kk,3);
            end
        end
    end
end

savevtk(V_shape, '/u/data/alcer/DFXRM_rec/Rec_test_2/V_weight.vtk');
savevtk(V_th, '/u/data/alcer/DFXRM_rec/Rec_test_2/V_th.vtk');
