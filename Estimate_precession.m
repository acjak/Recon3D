% This script locates, in the images summed at the different projections,
% the maximum intensity value. Aim: study misalignement between G vector
% and rotation axis

clear; close all;

% Load images from Fabio_isolate_blob
A = load('/u/data/alcer/DFXRM_rec/Rec_test_2/Sample2_cleaned.mat');

Data = A.foo;

% For each projection, find the position of the peak intensity
pos_max = zeros(size(Data,1), 3);
for i = 1:size(Data,1)
    if i < 98 | i > 121
        Layer = zeros(size(Data,2), size(Data,3));
        for j = 1:size(Data,2)
            for k = 1:size(Data,3)
                Layer(j,k) = Data(i,j,k);
            end
        end
    end
    s_l = smoothn(Layer);
    [I,J] = find(s_l == max(s_l(:)));
    pos_max(i,1) = i;
    pos_max(i,2) = I;
    pos_max(i,3) = J;
end

% Find the center of gravity of the points
CM_y = sum(pos_max(:,2))/size(pos_max,1);
CM_x = sum(pos_max(:,3))/size(pos_max,1);

% Central axis
P1 = [CM_x, CM_y, 0];
P2 = [CM_x, CM_y, 300];
pts = [P1; P2];

% Plot the motion of the intensity peak - 2D
figure; scatter(pos_max(:,3), pos_max(:,2), 'filled');
title('Intensity peak motion as the sample rotates');
xlabel('X (pixels)');
ylabel('Y (pixels)');
hold on;
% Plot CM
scatter(CM_x, CM_y, 'r', 'filled');
% Plot circle centered on the CM
r = 53;      % circle radius
c = [CM_x CM_y];  % circle center
n = 1000;   % number of points
%// running variable
t = linspace(0,2*pi,n);
x = c(1) + r*sin(t);
y = c(2) + r*cos(t);
%// draw line
line(x,y);
axis([50 200 120 240]);
axis equal

% Intensity peak in 3D
figure; scatter3(pos_max(:,3), pos_max(:,2), pos_max(:,1), 'filled');
title('Intensity peak motion as the sample rotates');
xlabel('X (pixels)');
ylabel('Y (pixels)');
zlabel('Projection number');
hold on;
scatter(CM_x, CM_y, 'r', 'filled');

% Intensity peak in 3D with rotation axis
figure; scatter3(pos_max(:,3), pos_max(:,2), pos_max(:,1), 'filled');
title('Intensity peak motion as the sample rotates');
xlabel('X (pixels)');
ylabel('Y (pixels)');
zlabel('projection number');
hold on;
plot3(pts(:,1), pts(:,2), pts(:,3), 'r', 'LineWidth',2);

% check coordinates
P1 = [CM_x, CM_y , 0];
P2 = [CM_x, CM_y, 15000];
pts = [P1; P2];
figure; mesh(squeeze(Data(1,:,:))); hold on; plot3(pts(:,1), pts(:,2), pts(:,3), 'r', 'LineWidth',2);
xlabel('X'); ylabel('Y');

fprintf('Estimated center of rotation: X = %03f, Y = %03f\n', CM_x, CM_y);
