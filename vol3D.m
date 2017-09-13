% This scipt compares the reconstructed 3D data with the sum of the
% intensities recorded at different projections

close all; clear;

addpath('/npy_matlab_master/');
% Read reconstructed volume. Format: X, Y, Z, param. Parameters: gamma, mu,
% completeness
Summed_img = readNPY('/u/data/alcer/DFXRM_rec/Rec_test_2/summed_data_astra.npy');
V = readNPY('/u/data/alcer/DFXRM_rec/Rec_test_2/grain_ang.npy');

%for i = [1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 130 135 140 145 150 155 160]
for i = 1:5:160
    proj_number = i;
    compare_shapes(proj_number, V, Summed_img);
end
    
function compare_shapes(alpha, V, Summed_img)

% Select volume corresponding to a certain completeness
V_th_low = zeros(size(V,1), size(V,2), size(V,1));
V_th_up = zeros(size(V,1), size(V,2), size(V,1));
for ii =1:size(V,1)
    for jj = 1:size(V,2)
        for kk = 1:size(V,3)
            if V(ii,jj,kk,3) > 0.5
                V_th_low(ii,jj,kk) = V(ii,jj,kk,3);
            end
            if V(ii,jj,kk,3) > 0.75
                V_th_up(ii,jj,kk) = V(ii,jj,kk,3);
            end
        end
    end
end

% Save calculated volumes
%savevtk(V_shape, '/u/data/alcer/DFXRM_rec/Rec_test_2/V_weight.vtk');
%savevtk(V_th, '/u/data/alcer/DFXRM_rec/Rec_test_2/V_th.vtk');

% Rotation angle
alpha1 = alpha*1.12;

C_tot = sum(sum(sum(V_th_low)));
% Calculate centre of mass of V_th_up, weighted by completeness
X_CM = 0; Y_CM = 0; Z_CM = 0;
for ii = 1:size(V,1)
    for jj = 1:size(V,2)
        for kk = 1:size(V,3)
            if V_th_low(ii,jj,kk) > 0.5
                X_CM = X_CM + ii*V_th_low(ii,jj,kk);
                Y_CM = Y_CM + jj*V_th_low(ii,jj,kk);
                Z_CM = Z_CM + kk*V_th_low(ii,jj,kk);
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
Sum_Z = sum(V_th_low,3);

% Plot binarized projection along Z, together with CM
%figure;
%Sum_Z(CM(1), CM(2)) = 100;
%h = pcolor(Sum_Z); shading flat;

% Find radius of the circle containing the binarized proj of the volume
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

Circle_container = zeros(int8(2*R) + 4, int8(2*R) + 4);
Dx = abs(int8(R) - CM(1));
Dy = abs(int8(R) - CM(2));
for ii = 3:size(V,1) + 2
    for jj = 3:size(V,2) + 2
        Circle_container(ii + Dx,jj + Dy) = Sum_Z(ii-2, jj-2);
    end
end

% Plot circle containing grain, and projected image before and after
% padding
figure;
subplot(1,2,1);
h = pcolor(Sum_Z); shading flat;
hold on;
th = 0:pi/50:2*pi;
xunit = R * cos(th) + CM(2);
yunit = R * sin(th) + CM(1);
plot(xunit, yunit);
hold on;
scatter(CM(2), CM(1), 'o', 'MarkerFaceColor', 'b');
title('Before padding');
subplot(1,2,2);
h = pcolor(Circle_container); shading flat;
hold on;
th = 0:pi/50:2*pi;
xunit = R * cos(th) + R + 3;
yunit = R * sin(th) + R + 3;
plot(xunit, yunit);
hold on;
scatter(CM(2) + Dy + 2, CM(1) + Dx + 2, 'o', 'MarkerFaceColor', 'b');
title('After padding');

IM_x_th_low = zeros(size(Circle_container,1), size(Circle_container,2)); 
IM_x_th_up = zeros(size(Circle_container,1), size(Circle_container,2));

Cyliner_container = zeros(size(Circle_container,1), size(Circle_container, 2), size(V,3));


% Sum reconstructed data along X (rotate and project the sample)
for ii = 1:size(Circle_container,1)-1
    for jj = 1:size(,2)-1
        for kk = 1:size(V,3)-1
            disp(ii * cosd(alpha) - jj * sind(alpha))
            disp(ii * sind(alpha) + jj * cosd(alpha))
            IM_x_th_low(jj,kk) = IM_x_th_low(jj,kk) + V_th_low((ii-CM(2)) * cosd(alpha) - (jj - CM(1)) * sind(alpha), (ii-CM(2)) * sind(alpha) + (jj-CM(1)) * cosd(alpha), kk);
            IM_x_th_up(jj,kk) = IM_x_th_up(jj,kk) + V_th_up((ii-CM(2)) * cosd(alpha) - (jj - CM(1)) * sind(alpha), (ii-CM(2)) * sind(alpha) + (jj-CM(1)) * cosd(alpha), kk);
        end
    end
end

IM_x_th_low = rot90(rot90(rot90(IM_x_th_low)));
IM_x_th_up = rot90(rot90(rot90(IM_x_th_up)));

% Find perimeter of the projected data and compare with experimental data
% (summed for a certain projection)
IM_x_bin_low = zeros(size(V,1)-1, size(V,2)-1);
IM_x_bin_up = zeros(size(V,1)-1, size(V,2)-1);

for ii = 1:size(V,1)-1
    for jj = 1:size(V,2)-1
        if IM_x_th_low(ii,jj) > 0
            IM_x_bin_low(ii,jj) = 1;
        end
        if IM_x_th_up(ii,jj) > 0
            IM_x_bin_up(ii,jj) = 1;
        end
    end
end

IM_x_th_low_r = imresize(IM_x_th_low, 3);
IM_x_bin_low_r = imresize(IM_x_bin_low, 3);

IM_x_th_up_r = imresize(IM_x_th_up, 3);
IM_x_bin_up_r = imresize(IM_x_bin_up, 3);

P1 = bwperim(IM_x_bin_low_r);
P3 = bwperim(IM_x_bin_up_r);

S1 = squeeze(Summed_img(:,alpha,:));
S2 = squeeze(Summed_img(:,80,:));

%figure; imagesc(S1);

overlay1 = flipud(imoverlay(imoverlay(S1, P1, [.3 1 .3]), P3, [.3 1 .3]));
%imshow(overlay1);

fig = figure; 
title(sprintf('Projection %03i', alpha));
subplot(1,3,1);
h = pcolor(S1); shading flat;
%title('Sum of diffraction images','fontsize',18);
subplot(1,3,2);
h = pcolor(IM_x_th_low_r); shading flat;
%title('Projected reconstruction','fontsize',18);
subplot(1,3,3);
imshow(overlay1); 
%title('Perimeter overlay','fontsize',18);
print(fig, sprintf('Shape_comp%03i', alpha), '-dpdf');
close;

%figure; imagesc(rot90(IM_x_shape));
%figure; imagesc(rot90(IM_y_shape));

end
