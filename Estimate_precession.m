% This script locates, in the images summed at the different projections,
% the centyer of mass of the intensity value. Aim: study horizontal shift
% of the rotation axis.

% Two cases are considered: CM using recorded intensity, and CM using
% position only (equivalent to consider a binary image)

clear; close all;

% Load images saved by img_sum.py
A = load('/u/data/alcer/DFXRM_rec/Rec_test_2/Sample2_cleaned.mat');

Data = A.foo;

% For each projection, find the position of the peak intensity
pos_max = zeros(size(Data,1), 3);
pos_max_bin = zeros(size(Data,1), 3);
for i = 1:size(Data,1)
    if i < 98 || i > 121
        Layer = zeros(size(Data,2), size(Data,3));
        for j = 1:size(Data,2)
            for k = 1:size(Data,3)
                Layer(j,k) = Data(i,j,k);
            end
        end
        SmL = smooth(Layer);
        X_CM = 0;
        X_CM_bin = 0;
        Y_CM = 0;
        Y_CM_bin = 0;
        for j = 1:size(Data,2)
            for k = 1:size(Data,3)
                if Layer(j,k) > 0
                    X_CM = X_CM + Layer(j,k) *j;
                    Y_CM = Y_CM + Layer(j,k)*k;
                    X_CM_bin = X_CM_bin + j; 
                    Y_CM_bin = Y_CM_bin + k; 
                    
                end
            end
        end
        pos_max(i,1) = i;
        pos_max(i,2) = X_CM/sum(sum(Layer));
        pos_max(i,3) = Y_CM/sum(sum(Layer));
        
        pos_max_bin(i,1) = i;
        pos_max_bin(i,2) = X_CM_bin/nnz(Layer);%sum(sum(Layer));
        pos_max_bin(i,3) = Y_CM_bin/nnz(Layer);%sum(sum(Layer));
    end
end

% Remove zero entries
pos_max( ~any(pos_max(:,2),2), : ) = [];
pos_max_bin( ~any(pos_max_bin(:,2),2), : ) = [];

% Find the center of gravity of the points
CM_y = sum(pos_max(:,2))/size(pos_max,1);
CM_x = sum(pos_max(:,3))/size(pos_max,1);
CM_y_bin = sum(pos_max_bin(:,2))/size(pos_max_bin,1);
CM_x_bin = sum(pos_max_bin(:,3))/size(pos_max_bin,1);

% Plot centre of mass for each projection, together with global CM
figure; 
subplot(1,2,1);
plot(pos_max(:,2), pos_max(:,3), '-.');
xlabel('X'); ylabel('Y');
hold on;
scatter(CM_y, CM_x, 500, '.r')
title('Weighted case');
subplot(1,2,2);
plot(pos_max_bin(:,2), pos_max_bin(:,3), '-.');
xlabel('X'); ylabel('Y');
hold on;
scatter(CM_y_bin, CM_x_bin, 500, '.r');
title('Binary case');
