clear all;
close all;
clc
A = imread('boat.bmp');
[im_ver, im_hor] = size(A);

% 以 DCT-II overcomplete dictionary
% 利用 OMP 做 sparse coding 來找出係數與相關的 atom

% 構造overcomplete DCT frame
N_1 = 8;  %每個Block 為8by8
N_2 = 16;  %Dictionary 由 16*16 個 8by8 2D-DCT 組成 
% DCT 1D
mDCT_oc = zeros(N_1, N_2);
for k = 1 : N_2
    if k == 1
        mDCT_oc(:, k) = (sqrt(1/N_1))*cos((pi/N_2)*(k-1)*((0:(N_1-1))' + (1/2)*ones(N_1,1)));
    else
        Temp = cos((pi/N_2)*(k-1)*((0:(N_1-1))' + (1/2)*ones(N_1,1)));
        Temp = Temp - mean(Temp);
        mDCT_oc(:, k) = Temp/norm(Temp);
    end
end

% 製造 DCT_2D overcomplete DCT dictionray
DCT_2D_frame = zeros(N_1*N_1, N_2*N_2);
for i = 1 : N_2 %vertical index
    for j = 1 : N_2 %horzontal index
        if j == 1
            Temp2D = mDCT_oc(:,i)*mDCT_oc(:,j)';
            DCT_2D_frame(:, (i-1)*N_2 + j) = Temp2D(:);
        else
        Temp2D = mDCT_oc(:,i)*mDCT_oc(:,j)';
        Temp_col = Temp2D(:) - mean(Temp2D(:))*ones(size(Temp2D(:))); %讓每個8by8的平均值為0
        DCT_2D_frame(:, (i-1)*N_2 + j) = Temp_col/norm(Temp_col); %normalizing
        end
    end
end
% DCT_2D overcomplete DCT dictionary 構造完畢

Sigma = 20;
A_n = double(A) + Sigma*randn(im_ver, im_hor); %產生含有雜訊的影像
Block_im_n = im2col(A_n, [8 8],'sliding'); %將含有雜訊的影像，用可以重疊的8by8 block 取出來，排成矩陣的column

%每個Block_im_n demean
Block_im_n_mean = mean(Block_im_n);
Block_im_n_mean0 = Block_im_n - kron(Block_im_n_mean, ones(64,1));


%% Denoising 程序開始!!
%階段1  用overcomplete DCT 求出每個分塊影像的 sparse表示係數
Block_sparse_coe = OMPdenoise(DCT_2D_frame, Block_im_n_mean0, Sigma); %這部要跑很久
 
%階段2  在已知每個含雜訊的分塊矩陣的表示係數的情況下，用close-form 求出denoise 影像 hat_A
%直接套用 close-form 公式
% close-form 看起來很複雜，事實上它只是把所有有重疊的block附近依照
lam = 0.1; %參數

%Denominator 為重疊的影像每個pixel 要除上的權重
%構造權重
Denominator = lam*ones(im_ver, im_hor);
for i = 1 : (im_ver-8+1)
    for j = 1 : (im_hor-8+1)
        Denominator(i:(i+7), j:(j+7)) = Denominator(i:(i+7), j:(j+7)) + ones(8,8);
    end
end
%構造有重疊的 sparse reconstruction
Recover_block_IM = DCT_2D_frame*Block_sparse_coe + Block_im_n_mean;
De_noise_im = lam*A_n;
coe_ind = 1; % Block_sparse_coe 的 index
for i = 1 : (im_ver-8+1)
    for j = 1 : (im_hor-8+1)
        De_noise_im(j:(j+7), i:(i+7)) = De_noise_im(j:(j+7), i:(i+7)) + reshape(Recover_block_IM(:, coe_ind), [8,8]);
        coe_ind = coe_ind + 1; 
    end
end

%% 秀出結果
figure;
subplot(1,3,1)
imshow(A);
xlabel(['PSNR = ', num2str(psnr(A,A))]);
title('Original Image');
subplot(1,3,2)
imshow(uint8(A_n));
xlabel(['PSNR = ', num2str(psnr(uint8(A_n),A))]);
title('Noising Image');
subplot(1,3,3);
Im_denoise = uint8(255*mat2gray(De_noise_im./Denominator));
imshow(Im_denoise);
xlabel(['PSNR = ', num2str(psnr(Im_denoise,A))])
title('Denoising Image');


