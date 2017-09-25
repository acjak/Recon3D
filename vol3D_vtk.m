% This script compares the reconstructed 3D data with the sum of the
% intensities recorded at different projections

close all; clear;

addpath('/npy_matlab_master/');
% Read reconstructed volume. Format: X, Y, Z, param. Parameters: gamma, mu,
% completeness
V = readNPY('/u/data/alcer/DFXRM_rec/Rec_test_2/grain_ang.npy');

V_th= zeros(size(V,1), size(V,2), size(V,1));
for ii =1:size(V,1)
    for jj = 1:size(V,2)
        for kk = 1:size(V,3)
            % The minimum completeness value for a voxel 
            % to be part of the volume is 0.5
            if V(ii,jj,kk,3) > 0.5
                V_th(ii,jj,kk) = V(ii,jj,kk,3);
            else
                V(ii,jj,kk,3) = NaN;
            end
        end
    end
end

% Save calculated volumes
savevtk(V_th, '/u/data/alcer/DFXRM_rec/Rec_test_2/V_th.vtk');
