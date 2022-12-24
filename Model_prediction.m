% This code performs model prediction for
% the human response to the perceived translucency 
% from the rating experiment.

clear all
close all

% the order of down sampling
smoothing = 1;
% range of excluded pixel from ROI
excluded_pix = 24;

% align the order of image features to that of human response data.
sort_ind =zeros(40,1);
cnt1=1;
cnt2=1;
for i = 1:40
    if rem(i,2) == 1
        sort_ind(cnt1) = i;
        cnt1 = cnt1 + 1;
    else
        sort_ind(cnt2+20) = i;
        cnt2 = cnt2 + 1;
    end
end
% load highlight and shading related features
features_highlight_shading = load('data/Features_1sd_ind.mat')
features_highlight_shading.Y_diff_mask = features_highlight_shading.Y_diff_mask(:,:,sort_ind);
features_highlight_shading.Y_spec_ind = features_highlight_shading.Y_spec_ind(:,:,sort_ind);
features_highlight_shading.diffuse_contrast = features_highlight_shading.diffuse_contrast(sort_ind);
features_highlight_shading.highlight_cov = features_highlight_shading.highlight_cov(sort_ind);
features_highlight_shading.mean_sharpness = features_highlight_shading.mean_sharpness(sort_ind);
[h, w, d] = size(features_highlight_shading.Y_diff_mask)

% load weights for the prediction model
load('data/linear_reg_n.mat')

% load luminance orientation map
dd_single_ori = load(['data/orientation_map_blurlvl_', num2str(smoothing), '.mat']);
dd_single_ori.ori_map_img = dd_single_ori.ori_map_img(sort_ind);

% extract highlights and shading related image features inside ROI  
for num = 1:40
    if num <=20
        features_highlight_shading.ori_map_img{num} = dd_single_ori.ori_map_img{num}(1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),:);
        features_highlight_shading.Y_spec_ind_cell{num} = features_highlight_shading.Y_spec_ind(1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),num)
        features_highlight_shading.Y_diff_mask_cell{num} = features_highlight_shading.Y_diff_mask(1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),num)
    else
        features_highlight_shading.ori_map_img{num}=dd_single_ori.ori_map_img{num}(1+1.5*excluded_pix*(smoothing+1):h-0.75*excluded_pix*(smoothing+1),1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),:);
        features_highlight_shading.Y_spec_ind_cell{num} = features_highlight_shading.Y_spec_ind(1+1.5*excluded_pix*(smoothing+1):h-0.75*excluded_pix*(smoothing+1),1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),num)
        features_highlight_shading.Y_diff_mask_cell{num} = features_highlight_shading.Y_diff_mask(1+1.5*excluded_pix*(smoothing+1):h-0.75*excluded_pix*(smoothing+1),1+excluded_pix*(smoothing+1):h-excluded_pix*(smoothing+1),num)
    end
end

% extract orientation anisotropy
imgnum = 1
while imgnum <= 40
    roinum_single_s = 1;
    roinum_single_d= 1;
    [xsize ysize dim] = size(features_highlight_shading.ori_map_img{imgnum});
    for i = 1:xsize
        for j = 1:ysize
            if features_highlight_shading.ori_map_img{imgnum}(i,j,2) ~= 0 && features_highlight_shading.Y_spec_ind_cell{imgnum}(i,j) == 1
                spec_aniso{imgnum}(roinum_single_s) = features_highlight_shading.ori_map_img{imgnum}(i,j,2);
                roinum_single_s = roinum_single_s +1;
            elseif features_highlight_shading.ori_map_img{imgnum}(i,j,2) ~= 0 && features_highlight_shading.Y_diff_mask_cell{imgnum}(i,j) ==1
                diff_aniso{imgnum}(roinum_single_d) = features_highlight_shading.ori_map_img{imgnum}(i,j,2);
                roinum_single_d = roinum_single_d +1;
            end
        end
    end
    imgnum = imgnum +1;
end

% calculate anisoshading ratio
for num = 1:40
    aniso_ratio(num) = mean(spec_aniso{num}(:)) / mean(diff_aniso{num}(:));
end

% z-score normalization
diffuse_contrast_all_n = zscore(features_highlight_shading.diffuse_contrast);
mag_grad_all_n = zscore(features_highlight_shading.mean_sharpness);
aniso_ratio_all_n = zscore(aniso_ratio);
highlight_cov_n = zscore(features_highlight_shading.highlight_cov);


% load the results of rating experiment
load('data/all_data_SSS_rating.mat');
all_res =All_data;

% model fitting
y_ss = fitall_n.Coefficients.Estimate(2) * diffuse_contrast_all_n + fitall_n.Coefficients.Estimate(3) * aniso_ratio_all_n + fitall_n.Coefficients.Estimate(4) * mag_grad_all_n + fitall_n.Coefficients.Estimate(5) * highlight_cov_n + fitall_n.Coefficients.Estimate(1);

% plot the result of prediction
figure;
scatter(zscore(y_ss),all_res ,'filled');
hold on
xlim([-3 3]);
ylim([-3 3]);
ax = gca;
ax.FontSize = 16;
h1 = plot([-2.5:2.5], [-2.5:2.5]);
hold off
h1.Color ='k'
xlabel('z-scored model output')
ylabel('z-scored rating value')

corr_coef = corrcoef(zscore(y_ss),all_res)