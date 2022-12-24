% This code calculates highlight sharpness, highlight coverage, and diffuse
% contrast.

clear all
close all

%% initial settings
% rendering parameter for my stimulus images
sss_num = [{'0.00'}, {'0.25'}, {'0.50'}, {'0.75'}, {'1.00'}];
freq_num = [{'10'}, {'20'}];
bump_num = [{'10'}, {'20'}];
angle_num = [{'45'}, {'80'}]

stimulus_cnt = 1;
stimulus_num = size(sss_num,2) * size(freq_num,2) * size(bump_num,2) * size(angle_num,2)
% define outer pixel range from ROI
pixel_range=48;
% coefficient for standard deviation to detect highlight regions
% If your stimulus image has a large standard deviation of luminance
% distribution, this value might be better set to larger than 1.
multiple_ind = 1;

mean_sharpness = zeros(1,stimulus_num);
highlight_cov= zeros(1,stimulus_num);
diffuse_contrast = zeros(1,stimulus_num);
for sss = 1:5
    for freq = 1:2
        for bump = 1:2
            for angle = 1:2
                
                % load image
                img = imread(['Images/SSS_', sss_num{sss}, '_freq_0.', freq_num{freq}, '_meso_0.',bump_num{bump},'_angle_-',angle_num{angle}, '.bmp']);
                img = imresize(img,1/2)
                [height, width, dim] = size(img);
                if stimulus_cnt == 1
                    Y_spec_ind = zeros(height,width,stimulus_num);
                    Y_diff_mask = zeros(height,width,stimulus_num);
                end
                % convert sRGB to xyY color space
                R_v = reshape(img(:,:,1), [1, height * width]);
                G_v = reshape(img(:,:,2), [1, height * width]);
                B_v = reshape(img(:,:,3), [1, height * width]);
                sRGB_v = [R_v;G_v;B_v];
                RGB_v = SRGBGammaUncorrect(sRGB_v);
                XYZ_v = SRGBPrimaryToXYZ(RGB_v);
                xyY_v = XYZToxyY(XYZ_v);
                xyY_v(1,:) = 0.313;
                xyY_v(2,:) = 0.329;
                
                % extract lumiannce information
                Y_data_v = xyY_v(3,:);
                Y_data = reshape(Y_data_v,[height, width]);
                
                % calculate mean and standard deviation of luminance
                % distribution on surface of the image
                mean_Y = mean(Y_data_v(Y_data_v ~= 0))
                std_Y = std(Y_data_v(Y_data_v ~= 0))
                highlight_ind = mean_Y + multiple_ind * std_Y ;
                
                % detect pixels where luminance edge
                [BW1, threshold] = edge(Y_data,'Canny', [0.02, 0.04]);
                BW1 = BW1(pixel_range:(height-(pixel_range+1)), pixel_range:(width-(pixel_range+1)));
                
                % smoothing image and calculate highlight sharpness
                [dx, dy] = smooth_gradient(Y_data, 0.25);
                dx = dx(pixel_range:(height-pixel_range+1), pixel_range:(width-pixel_range+1));
                dy = dy(pixel_range:(height-pixel_range+1), pixel_range:(width-pixel_range+1));
                sharpness(:,:,stimulus_cnt) = hypot(dx, dy);
                
                [target_row, target_col]= find(BW1==1);
                target_sharpness = zeros(size(target_row,1),1)
                for i = 1:size(target_row,1)
                    target_sharpness(i) = sharpness(target_row(i),target_col(i),stimulus_cnt);
                end
                mean_sharpness(stimulus_cnt) = mean(target_sharpness(:));
                
                % devide the surface into the specular highlight regions and the diffuse
                % shading regions
                for i = 1:width
                    for j = 1:height
                        
                        if Y_data(i,j) < highlight_ind && Y_data(i,j) > 0.001
                            Y_diff_mask(i,j,stimulus_cnt) = 1;
                        elseif Y_data(i,j) >= highlight_ind && Y_data(i,j) > 0.001
                            Y_spec_ind(i,j, stimulus_cnt) = 1;
                        end
                    end
                end
                
                % calculate highlight coverage
                highlight_cov(stimulus_cnt) = size(find(Y_spec_ind(:,:,stimulus_cnt) == 1),1)  
                
                % calculate diffuse contrast
                Y_diff_mask_v = reshape(Y_diff_mask(:,:,stimulus_cnt),[1,width*height]);
                target_Y = Y_data_v( Y_diff_mask_v == 1);
                diffuse_contrast(stimulus_cnt) = std(target_Y)/mean(target_Y);
                stimulus_cnt = stimulus_cnt + 1;
            end
        end
    end
end

save('data/Features_1sd_ind.mat', 'diffuse_contrast', 'highlight_cov', 'mean_sharpness','Y_spec_ind', 'Y_diff_mask')

%% define smooth_gradient
function [GX, GY] = smooth_gradient(I, sigma)

filterExtent = ceil(4*sigma);
x = -filterExtent:filterExtent;

c = 1/(sqrt(2*pi)*sigma);
gaussKernel = c * exp(-(x.^2)/(2*sigma^2));
gaussKernel = gaussKernel/sum(gaussKernel);
GaussKernel = gradient(gaussKernel);
negVals = GaussKernel < 0;
posVals = GaussKernel > 0;
GaussKernel(posVals) = GaussKernel(posVals)/sum(GaussKernel(posVals));
GaussKernel(negVals) = GaussKernel(negVals)/abs(sum(GaussKernel(negVals)));

GX = imfilter(I, GaussKernel, 'conv', 'replicate');
GY  = imfilter(I, GaussKernel', 'conv', 'replicate');
end