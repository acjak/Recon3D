function blobs_layer = Blobs_finder(npy_array, mean_value, z_value, idx_tgw)
% This function returns, for each Z value of an npy array, the regions that
% should be associated with the shape of a grain
% npy_array from recon3d.py
% mean_value is the max(mean) of the weight, calculated for each layer
% z_value is considered layer
% idx_wtg is the voxel index to consider. 1: theta, 2: gamma, 3: weight

% Play with a single slice
if idx_tgw == 3
    B = (npy_array(:, :, z_value, idx_tgw)/mean_value);
elseif idx_tgw == 1
    B = npy_array(:, :, z_value, idx_tgw);
end
D = (imregionalmin(B,8));     % Find local minima
Min_matrix = zeros(size(B));

for i = 1:size(B,1)
    for j = 1:size(B,2)
        if D(i,j) == 1
            Min_matrix(i,j) = B(i,j);
        end
    end
end

% Higher local minima
D_max = max(max(Min_matrix));
% Lower local minima
Min_matrix(~Min_matrix) = inf;
D_min = min(min(Min_matrix));

% Find the area of the basins relative to the minima
n_steps = 1000;  % Number of steps between min and max
prop_steps = zeros(1,4);
for aa = D_min:(1/n_steps):D_max
    B_one_level = zeros(size(B));
    for bb = 1:size(B,1)
        for cc = 1:size(B,2)
            if B(bb,cc) < aa
                B_one_level(bb,cc) = 1;
            end
        end
    end
    % Find the different regions in the image
    L2 = labelmatrix(bwconncomp(B_one_level));
    n_reg = size(unique(L2),1) - 1;
    for dd = 1:n_reg
        L3 = zeros(size(L2));
        for ee = 1:size(L2,1)
            for ff = 1:size(L2,2)
                if L2(ee,ff) == dd
                    L3(ee,ff) = 1;
                end
            end
        end
        % Calculate area and centroid of the considered region
        s = regionprops(L3, 'Centroid', 'Area');
        if s.Area > 50 && s.Area < 400
            % Use a frame (the grain is well centred)
            if s.Centroid(1) > 20 && s.Centroid(1) < 80
                if s.Centroid(2) > 20 && s.Centroid(2) < 80
                    prop_steps = [prop_steps; [aa s.Centroid(1) s.Centroid(2) s.Area]];
                end
            end
        end
    end
end

prop_steps( ~any(prop_steps,2), : ) = [];

% Look at the CM distribution and, for each cluster, find minima corresponding
% to largest area
Loc_CM = zeros(size(B));
for u = 1:size(prop_steps,1)
    Loc_CM(int8(prop_steps(u, 2)), int8(prop_steps(u,3))) = 1;
end
%figure; imshow(Loc_CM);
Loc_CM_lab = bwlabel(Loc_CM); % label the different clusters of minima

n1 = size(unique(Loc_CM_lab),1) - 1;
val_max = zeros(n1, 5);
for uu = 1:n1
    M1 = zeros(size(B));
    for i = 1:size(B,1)
        for j = 1:size(B,2)
            if Loc_CM_lab(i,j) == uu
                M1(i,j) = 1;
            end
        end
    end

    Data_M1 = zeros(1,4);
    for vv = 1:size(prop_steps, 1)
        if M1(round(prop_steps(vv,2)), round(prop_steps(vv,3))) == 1
            Data_M1 = [Data_M1; prop_steps(vv,:)];
        end
    end
    Data_M1( ~any(Data_M1,2), : ) = [];
    [max_A, idx_max_A] = max(Data_M1(:,4));
    val_max(uu,1:4) = Data_M1(idx_max_A,:); % aa x_min y_min A_ min
    val_max(uu,5) = idx_max_A;     % min position in prop_steps
end

% Save, in a single image, the blobs corresponding to the local minima

Blobs_final_comb = zeros(size(B));
for uu = 1:n1
    for aa = val_max(uu,1)
        B_one_level = zeros(size(B));
        for bb = 1:size(B,1)
            for cc = 1:size(B,2)
                if B(bb,cc) < aa
                    B_one_level(bb,cc) = 1;
                end
            end
        end
        % Find the different regions in the image
        L2 = labelmatrix(bwconncomp(B_one_level));
        n_reg = size(unique(L2),1) - 1;
        for dd = 1:n_reg
            L3 = zeros(size(L2));
            for ee = 1:size(L2,1)
                for ff = 1:size(L2,2)
                    if L2(ee,ff) == dd
                        L3(ee,ff) = 1;
                    end
                end
            end
            % Calculate area and centroid of the considered region
            s = regionprops(L3, 'Centroid', 'Area');
            if s.Centroid(1) == val_max(uu,2) && s.Centroid(2) == val_max(uu,3) && s.Area == val_max(uu,4)
                Blobs_final_comb = Blobs_final_comb + L3;
            end
        end
    end
end

blobs_layer = Blobs_final_comb;

end  % function
