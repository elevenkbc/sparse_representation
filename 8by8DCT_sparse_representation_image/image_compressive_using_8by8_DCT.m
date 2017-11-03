clear all;
close all;
clc
A = imread('barbara.bmp');
[im_ver, im_hor] = size(A);
block_ver_num = im_ver/8;
block_hor_num = im_hor/8;

% mDCT8 = zeros(8, 8); %每一個column 為DCT-IV 的在 frequency space 的 basis
% for k = 1 : 8
%     mDCT8(:, k) = sqrt(2/8)*cos((pi/8)*((0:7)' + (1/2)*ones(8,1))*(k - 1/2));
% end

mDCT8 = zeros(8, 8); %每一個column 為DCT 在 frequency space 的 basis
for k = 1 : 8
    if k == 1
        mDCT8(:, k) = (1/2)*(1/sqrt(2))*cos((k-1)*pi*((2*(0:7)'+ones(8,1))/16));
    else
        mDCT8(:, k) = (1/2)*cos((k-1)*pi*((2*(0:7)'+ones(8,1))/16));
    end
    
end

%plot 8by8 2D DCT basis

DCTbasis = zeros(64,64);
for i = 1 : 8 %vertical index
    for j = 1 : 8 %horzontal index
        block_ind_start_i = (i-1)*8 + 1;
        block_ind_start_j = (j-1)*8 + 1;
        DCTbasis((block_ind_start_i:block_ind_start_i+7), (block_ind_start_j:block_ind_start_j+7)) = mDCT8(:,i)*mDCT8(:,j)';
    end
end
plot_8by8_grid(DCTbasis);
title('DCT  2D basis')

DCT_coe = zeros(im_ver, im_hor);
compressive_im = zeros(im_ver, im_hor);

%對於每一個 8by8 Block 做2D DCT
doubleA = double(A); %轉成 double 才能計算 DCT
for i = 1 : block_ver_num
    for j = 1 : block_hor_num
        %提出8by8 block
        block_ind_start_i = (i-1)*8 + 1;
        block_ind_start_j = (j-1)*8 + 1;
        temp_block = doubleA((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7));
        
        temp_dct = (mDCT8')*temp_block*mDCT8; %做 8by8 2D DCT
        DCT_coe((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7)) = temp_dct; %將係數塞回去
        % DCT 係數 threadhold
%         thrn = 7;
%         Thr_matrix = fliplr(toeplitz(zeros(1,8), [zeros(1, 8-thrn), ones(1,thrn)]));
%         
%         temp_dct_thr = temp_dct.*Thr_matrix;
        
        temp_dct_thr = temp_dct.*(abs(temp_dct) > 15);
        %IDCT
        
        compressive_im((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7)) = mDCT8*temp_dct_thr*mDCT8';
    end
end

figure('position', [100, 400, 1100, 350])
subplot(1,3,1)
imshow(mat2gray(DCT_coe));
title('DCT 係數')

subplot(1,3,2); imshow(A); colormap('gray')
title('原始圖片')
Ps = psnr(uint8(compressive_im), A);

subplot(1,3,3); imshow(uint8(compressive_im)); colormap('gray')
title('compressive 圖片')
xlabel(['PSNR = ', num2str(Ps)]);



