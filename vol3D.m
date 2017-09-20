% This scipt compares the reconstructed 3D data with the sum of the
% intensities recorded at different projections

close all; clear;

addpath('/npy_matlab_master/');
% Read reconstructed volume. Format: X, Y, Z, param. Parameters: gamma, mu,
% completeness
Summed_img = readNPY('/u/data/alcer/DFXRM_rec/Rec_test_2/summed_data_astra.npy');
V = readNPY('/u/data/alcer/DFXRM_rec/Rec_test_2/grain_ang.npy');

for i = [1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 130 135 140 145 150 155 160]
    proj_number = i;
    disp(i);
    [P_05, P_07, P_09] = compare_shapes(proj_number, V);

    P_05 = flipud(rot90(imresize(P_05(1:100, 1:100),3)));
    P_07 = flipud(rot90(imresize(P_07(1:100, 1:100),3)));
    P_09 = flipud(rot90(imresize(P_09(1:100, 1:100),3)));
    
    P_05_bin = zeros(size(P_05));
    P_07_bin = zeros(size(P_07));
    P_09_bin = zeros(size(P_09));
    
    for ii = 1:size(P_05,1)
        for jj = 1:size(P_05,2)
            if P_05(ii,jj) > 0
                P_05_bin(ii,jj) = 1;
            end
            if P_07(ii,jj) > 0
                P_07_bin(ii,jj) = 1;
            end
            if P_09(ii,jj) > 0
                P_09_bin(ii,jj) = 1;
            end
        end
    end
    
    % Fill holes in the binary images
    P_05_bin = imfill(P_05_bin);
    P_07_bin = imfill(P_07_bin);
    P_09_bin = imfill(P_09_bin);
    
    % Find perimeter of binary images
    Per_05 = bwperim(P_05_bin);
    Per_07 = bwperim(P_07_bin);
    Per_09 = bwperim(P_09_bin);
    
    % Select the summed intensities image corresponding to the considered
    % angle
    Sum = squeeze(Summed_img(:,proj_number,:));
    
    overlay1 = flipud(imoverlay(imoverlay(imoverlay(Sum, Per_05, [.3 1 .3]), Per_07, [.3 1 .3]), Per_09, [.3 1 .3]));
    
    F = figure; 
    subplot(1,3,1);
    h = pcolor(Sum); shading flat; hold on;
    subplot(1,3,2);
    h = pcolor(P_05); shading flat; hold on;
    subplot(1,3,3);
    imshow(overlay1);
    address_F = sprintf('Shape_comp/Shape_comp%03i.png', proj_number);
    saveas(F, address_F, 'png');
    close;
end
    
% For a given completeness value, and a given rotation angle, this function 
% calculates the shape of the diffraction signal collected by the detector, 
% and the border of the shape
function [Proj_final_05, Proj_final_07, Proj_final_09] = compare_shapes(alpha, V)

% Select volume corresponding to the lower completeness value considered (0.5)
V_th= zeros(size(V,1), size(V,2), size(V,1));
for ii =1:size(V,1)
    for jj = 1:size(V,2)
        for kk = 1:size(V,3)
            if V(ii,jj,kk,3) > 0.5
                V_th(ii,jj,kk) = V(ii,jj,kk,3);
            end
        end
    end
end

% Save calculated volumes
%savevtk(V_th, '/u/data/alcer/DFXRM_rec/Rec_test_2/V_th.vtk');

% Rotation angle
alpha1 = alpha*1.125;

C_tot = sum(sum(sum(V_th)));
% Calculate centre of mass of V_th_up, weighted by completeness
X_CM = 0; Y_CM = 0; Z_CM = 0;
for ii = 1:size(V,1)
    for jj = 1:size(V,2)
        for kk = 1:size(V,3)
            % Set a threshold to accept a voxel
            if V_th(ii,jj,kk) > 0.5
                X_CM = X_CM + ii*V_th(ii,jj,kk);
                Y_CM = Y_CM + jj*V_th(ii,jj,kk);
                Z_CM = Z_CM + kk*V_th(ii,jj,kk);
            end
        end
    end
end
CM = [X_CM, Y_CM, Z_CM] / C_tot;

% Define rotating sample and project it along X
% Start by determining a cylinder that contains the sample. First, project
% sample along Z
Sum_Z = zeros(size(V,1), size(V,2));
% Project all values along Z
Sum_Z = sum(V_th,3);

% Plot binarized projection along Z, together with CM
%figure;
Sum_Z(int8(CM(1)), int8(CM(2))) = 100;

% Find radius of the circle containing the entire XY projection
count_dist = zeros(nnz(Sum_Z), 2);
cc = 0;
for ii = 1:size(V,1)
    for jj = 1:size(V,2)
        if Sum_Z(ii,jj) > 0
            cc = cc + 1;
            count_dist(cc,1) = cc;
            count_dist(cc,2) = sqrt((ii-CM(1)).^2 + (jj-CM(2)).^2);
        end
    end
end
% Find radius containg all values
R = max(count_dist(:,2));

% Translate the projection, so that the CM is at the image center 
Circle_container = zeros(int8(2*R) + 4, int8(2*R) + 4);
Dx = abs(int8((size(Circle_container,1)/2) - CM(1)));
Dy = abs(int8((size(Circle_container,2)/2) - CM(2)));
Dz = abs(int8(size(V,3)/2 - CM(3)));
for ii = 1:size(V,1)
    for jj = 1:size(V,2)
        Circle_container(ii + Dx, jj + Dy) = Sum_Z(ii, jj);
    end
end

X_CM_circ = 0; Y_CM_circ = 0;
for aa = 1:size(Circle_container, 1)
    for bb = 1:size(Circle_container, 2)
        X_CM_circ = X_CM_circ + aa*Circle_container(aa,bb);
        Y_CM_circ = Y_CM_circ + bb*Circle_container(aa,bb);
    end
end
C_tot_circ = sum(sum(sum(V_th)));
CM_circ = [X_CM_circ, Y_CM_circ] / C_tot_circ;

% Plot projection before and after translation
% figure; 
% subplot(1,2,1);
% h = pcolor(Sum_Z); shading flat; hold on;
% scatter(CM(2), CM(1));
% subplot(1,2,2);
% h = pcolor(Circle_container); shading flat; hold on;
% scatter(CM_circ(2), CM_circ(1));

IM_x_th_low = zeros(size(Circle_container,1), size(Circle_container,2)); 
IM_x_th_up = zeros(size(Circle_container,1), size(Circle_container,2));

% Iterate for each layer. Z translation to make CM at the center of the
% projected image
Cylinder_container = zeros(size(Circle_container,1), size(Circle_container, 2), size(V,3) + Dz);
for ii = 3:size(V,1) + 2
    for jj = 3:size(V,2) + 2
        for kk = 3:size(V,3) + 2
            Cylinder_container(ii + Dx,jj + Dy, kk + Dz) = V_th(ii-2, jj-2, kk - 2);
        end
    end
end

% Rotate each horizontal layer by alpha and project (rotation around Z)
Rot_Z = zeros(size(imrotate(Cylinder_container, alpha1), 1), ...
    size(imrotate(Cylinder_container, alpha), 1), size(V,3));
Rot_Z = imrotate(Cylinder_container, alpha1+180);
 
% Rotate the sample around Y
theta = 10.3754;        % Scattering angle
test_l = squeeze(Rot_Z(:, 50, :));
test_l_rot = imrotate(test_l, theta); 

Rot_Y = zeros(size(test_l_rot, 1), size(Rot_Z,1), size(test_l_rot,2));
for jj = 1:size(Rot_Y,2)
    layer_in = squeeze(Rot_Z(:,jj,:));
    layer_fin = imrotate(layer_in, 10.3754);    
    for ii = 1:size(Rot_Y,1)
        for kk = 1:size(Rot_Y,3)
            Rot_Y(ii,jj,kk) = layer_fin(ii,kk);
        end
    end
end

% Find CM in Rot_Y (to check that it's approx. at the center of the volume)
X_CM_rot = 0; Y_CM_rot = 0; Z_CM_rot = 0;
C_tot_rot = sum(sum(sum(Rot_Y)));
for ii = 1:size(Rot_Y,1)
    for jj = 1:size(Rot_Y,2)
        for kk = 1:size(Rot_Y,3)
            if Rot_Y(ii,jj,kk) > 0
                X_CM_rot = X_CM_rot + ii*Rot_Y(ii,jj,kk);
                Y_CM_rot = Y_CM_rot + jj*Rot_Y(ii,jj,kk);
                Z_CM_rot = Z_CM_rot + kk*Rot_Y(ii,jj,kk);
            end
        end
    end
end
CM_rot = [X_CM_rot, Y_CM_rot, Z_CM_rot] / C_tot_rot;
     
% Sum the rotated volume along the X axis
X_sum_05 = zeros(size(Rot_Y,2), size(Rot_Y,3));
X_sum_05 = squeeze(sum(Rot_Y, 1));

% Resize the image, so that it matches the original size
Proj_final_05 = zeros(size(Sum_Z));
frame_x = int8((size(Rot_Y,2) - size(Proj_final_05,1))/2);
frame_y = int8((size(Rot_Y,3) - size(Proj_final_05,2))/2);
for aa = 1:size(Proj_final_05, 1)
    for bb = 1:size(Proj_final_05, 2)
        Proj_final_05(aa, bb) = X_sum_05(aa + frame_x, bb + frame_y);  
    end
end

% We also want to project the volumes corresponding to C = 0.7 and C = 0.9
Rot_Y_07 = zeros(size(Rot_Y));
Rot_Y_09 = zeros(size(Rot_Y));
for ii = 1:size(Rot_Y,1)
    for jj = 1:size(Rot_Y,2)
        for kk = 1:size(Rot_Y,3)
            if Rot_Y(ii,jj,kk) > 0.7
                Rot_Y_07(ii,jj,kk) = Rot_Y(ii,jj,kk);
            end
            if Rot_Y(ii,jj,kk) > 0.9
                Rot_Y_09(ii,jj,kk) = Rot_Y(ii,jj,kk);
            end
        end
    end
end

% Sum the rotated volume along the X axis
X_sum_07 = zeros(size(Rot_Y,2), size(Rot_Y,3));
X_sum_07 = squeeze(sum(Rot_Y_07, 1));

X_sum_09 = zeros(size(Rot_Y,2), size(Rot_Y,3));
X_sum_09 = squeeze(sum(Rot_Y_09, 1));

% Resize the image, so that it matches the original size
Proj_final_07 = zeros(size(Sum_Z));
frame_x = int8((size(Rot_Y,2) - size(Proj_final_07,1))/2);
frame_y = int8((size(Rot_Y,3) - size(Proj_final_07,2))/2);
for aa = 1:size(Proj_final_07, 1)
    for bb = 1:size(Proj_final_07, 2)
        Proj_final_07(aa, bb) = X_sum_07(aa + frame_x, bb + frame_y);  
    end
end

Proj_final_09 = zeros(size(Sum_Z));
frame_x = int8((size(Rot_Y,2) - size(Proj_final_09,1))/2);
frame_y = int8((size(Rot_Y,3) - size(Proj_final_09,2))/2);
for aa = 1:size(Proj_final_09, 1)
    for bb = 1:size(Proj_final_09, 2)
        Proj_final_09(aa, bb) = X_sum_09(aa + frame_x, bb + frame_y);  
    end
end

end
