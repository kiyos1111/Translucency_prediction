% This code calculates orientation field of luminance distribution on the stimulus images.
% Output is a 3-dimensional matrix.
% The first dim contains the direction of which the filter response is maximum.
% The second dim contains local orientation anisotropy of luminance distributions.
% The third dim contains binary info showing object region or background region. 

clear all
close all

%% initial settings
% rendering parameter for my stimulus images.
sss_num = [{'0.00'}, {'0.25'}, {'0.50'}, {'0.75'}, {'1.00'}];
freq_num = [{'10'}, {'20'}];
bump_num = [{'10'}, {'20'}];
angle_num = [{'45'}, {'80'}];

% load orientation filters
load('data/orientation_filters.mat');
OFs = extracted_filters;
blurlevel = 1;

imgnum=1;
for sss = 1:5
    for freq = 1:2
        for bump = 1:2
            for angle = 1:2
                img = imread(['Images/SSS_', sss_num{sss}, '_freq_0.', freq_num{freq}, '_meso_0.',bump_num{bump},'_angle_-',angle_num{angle}, '.bmp']);
                
                
                % calculate luminance values 
                [xsize ysize] = size(img(:,:,1));
                img_v = reshape(img,[xsize*ysize,3]);
                RGBimg_v = SRGBGammaUncorrect(img_v');
                XYZimg_v = SRGBPrimaryToXYZ(RGBimg_v);
                Ydata_v = XYZimg_v(2,:);
                Ydata = reshape(Ydata_v,[xsize, ysize]);
                Ydata(find(Ydata<0.001)) = 0;
                Y = Ydata; % This is calculated luminance values.
                
                % apply orientation filters to the stimulus images
                for i = 1:24
                    ori_filtered_Y(:,:,i) = imfilter(Y,OFs(:,:,i));
                end
                
                [xsize ysize] = size(ori_filtered_Y(:,:,1))
                ori_map_img{imgnum} = zeros(xsize/(2^blurlevel),ysize/(2^blurlevel));
                ori_filt_Y_sq = ori_filtered_Y.^2;
                
                for i = 1:24
                    filtered_Y(:,:,i) = blurDn(ori_filt_Y_sq(:,:,i), blurlevel);
                                    end
                mask = blurDn(Y,blurlevel);
                
                % detect maximum response filter
                for i = 1:xsize/(2^blurlevel)
                    for j = 1:ysize/(2^blurlevel)
                        if mask(i,j) ~= 0
                            ori_ind(i,j) = find(filtered_Y(i,j,:) == max(filtered_Y(i,j,:)),1);
                            ori_map_img{imgnum}(i,j,1) = ori_ind(i,j)./24*360;
                        else
                            ori_ind(i,j) = 0;
                            ori_map_img{imgnum}(i,j,1:3) = 0;
                        end
                        
                    end
                end
                
                cnt = 1;
                % calculate local luminance anisotropy
                for i = 1:xsize/(2^blurlevel)
                    for j = 1:ysize/(2^blurlevel)
                        if mask(i,j) ~= 0
                            A(cnt)=1-sqrt(min(filtered_Y(i,j,:)) / max(filtered_Y(i,j,:)));
                            ori_map_img{imgnum}(i,j,2) = A(cnt);
                            ori_map_img{imgnum}(i,j,3) = 1;
                            cnt = cnt +1;
                        end
                        
                    end
                end
   
                imgnum = imgnum+1;
            end
        end
    end
end

save(['data/orientation_map_blurlvl_', num2str(blurlevel), '.mat'],'ori_map_img');
