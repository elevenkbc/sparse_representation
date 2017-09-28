function plot_8by8_grid(im)
im = double(im);
[im_ver, im_hor] = size(im);
block_ver_num = im_ver/8;
block_hor_num = im_hor/8;

%�N�v��gray �� RGB  �åB�վ�� 0 ~ 255����
im_rgb = zeros(im_ver, im_hor, 3);
im_rgb(:,:,1) = mat2gray(im)*255;
im_rgb(:,:,2) = mat2gray(im)*255;
im_rgb(:,:,3) = mat2gray(im)*255;



%�N�v���X�j
im_a = zeros(im_ver+block_ver_num-1, im_hor+block_hor_num-1, 3);
%�N�X�j���v����l�Ƭ�����
im_a(:,:,1) = 255;


%�N�v����J�X�j���v��
for i = 1 : block_ver_num
    for j = 1 : block_ver_num
        block_indi_start = (i-1)*9+1;
        block_indj_start = (j-1)*9+1;
        
        im_a((block_indi_start:block_indi_start+7), (block_indj_start:block_indj_start+7), :) = im_rgb((((i-1)*8+1):((i-1)*8+8)),(((j-1)*8+1):((j-1)*8+8)), :);
    end
end

figure;
imshow(uint8(im_a));

end
