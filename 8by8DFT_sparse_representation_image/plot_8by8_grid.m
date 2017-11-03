function plot_8by8_grid(im)
im = double(im);
[im_ver, im_hor] = size(im);
block_ver_num = im_ver/8;
block_hor_num = im_hor/8;

%將影像gray 轉 RGB  並且調整到 0 ~ 255之間
im_rgb = zeros(im_ver, im_hor, 3);
im_rgb(:,:,1) = mat2gray(im)*255;
im_rgb(:,:,2) = mat2gray(im)*255;
im_rgb(:,:,3) = mat2gray(im)*255;



%將影像擴大
im_a = zeros(im_ver+block_ver_num+1, im_hor+block_hor_num+1, 3);
%將擴大的影像初始化為紅色
im_a(:,:,3) = 255;


%將影像填入擴大的影像
for i = 1 : block_ver_num
    for j = 1 : block_ver_num
        block_indi_start = (i-1)*9+2;
        block_indj_start = (j-1)*9+2;
        
        im_a((block_indi_start:block_indi_start+7), (block_indj_start:block_indj_start+7), :) = im_rgb((((i-1)*8+1):((i-1)*8+8)),(((j-1)*8+1):((j-1)*8+8)), :);
    end
end

figure;
imshow(uint8(im_a));
% hold on
% line_length_ver = size(im_a, 1) + 0.5;
% line_length_hor = size(im_a, 2) + 0.5;
% %畫直線
% for i = 1 : block_ver_num-1
%     plot([9*i, 9*i], [0.5, line_length_ver], 'r', 'linewidth', 2.5);
% end
% %畫橫線
% for i = 1 : block_hor_num-1
%     plot([0.5, line_length_hor], [9*i, 9*i], 'r', 'linewidth', 2.5);
% end
% hold off

end