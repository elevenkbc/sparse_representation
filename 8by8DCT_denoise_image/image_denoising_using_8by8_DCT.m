clear all;
close all;
clc
A = imread('barbara.bmp');
A_n = imread('barbara_n.bmp');
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

plot_8by8_grid(DCTbasis)

DCT_coe = zeros(im_ver, im_hor);
compressive_im = zeros(im_ver, im_hor);

%對於每一個 8by8 Block 做2D DCT
doubleA = double(A); %轉成 double 才能計算 DCT

noise_imag = double(A_n);
for i = 1 : block_ver_num
    for j = 1 : block_hor_num
        %提出8by8 block
        block_ind_start_i = (i-1)*8 + 1;
        block_ind_start_j = (j-1)*8 + 1;
        temp_block = noise_imag((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7));
        
        temp_dct = (mDCT8')*temp_block*mDCT8; %做 8by8 2D DCT
        DCT_coe((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7)) = temp_dct; %將係數塞回去
        % DCT 係數 threadhold
        temp_dct_thr = temp_dct.*(abs(temp_dct) > 26);
        %IDCT
        
        compressive_im((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7)) = mDCT8*temp_dct_thr*mDCT8';
    end
end

figure('Position',[100,300,1100,350]) 
subplot(1,3,1); imshow(A); title('原始圖片');
subplot(1,3,2); imshow(A_n); title('含有雜訊圖片');
Ps1 = psnr(A_n, A);
xlabel(['PSNR = ', num2str(Ps1)]);
subplot(1,3,3); imshow(uint8(compressive_im)); title('denoise by DCT');
Ps2 = psnr(uint8(compressive_im), A);
xlabel(['PSNR = ', num2str(Ps2)]);



