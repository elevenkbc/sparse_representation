function image_denoising_by_assign_oc_frame(assign_frame, original_im, noising_im, Sigma)
d_noising_im = double(noising_im);
[im_ver, im_hor] = size(d_noising_im);
Block_im_n = im2col(d_noising_im, [8, 8],'sliding'); %將含有雜訊的影像，用可以重疊的8by8 block 取出來，排成矩陣的column
%每個Block_im_n 去除平均值
%並且將每個 Block 的平均值存成 Block_im_n_mean
Block_im_n_mean = mean(Block_im_n);
Block_im_n_mean0 = Block_im_n - kron(Block_im_n_mean, ones(64,1));


%% Denoising 程序開始!!
%階段1  用overcomplete DCT 求出每個分塊影像的 sparse表示係數
Block_sparse_coe = OMPdenoise(assign_frame, Block_im_n_mean0, Sigma); %這部要跑很久
 
%階段2  在已知每個含雜訊的分塊矩陣的表示係數的情況下，用close-form 求出denoise 影像 hat_A
%直接套用 close-form 公式
% close-form 看起來很複雜，事實上它只是把所有有重疊的block附近依照
lam = 0.3; %參數

%Denominator 為重疊的影像每個pixel 要除上的權重
%構造權重
Denominator = lam*ones(im_ver, im_hor);
for i = 1 : (im_ver-8+1)
    for j = 1 : (im_hor-8+1)
        Denominator(i:(i+7), j:(j+7)) = Denominator(i:(i+7), j:(j+7)) + ones(8,8);
    end
end
%構造有重疊的 sparse reconstruction
Recover_block_IM = assign_frame*Block_sparse_coe + Block_im_n_mean;
De_noise_im = lam*d_noising_im;
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
imshow(original_im);
xlabel(['PSNR = ', num2str(psnr(original_im,original_im))]);
title('Original Image');
subplot(1,3,2)
imshow(uint8(d_noising_im));
xlabel(['PSNR = ', num2str(psnr(uint8(d_noising_im),original_im))]);
title('Noising Image');
subplot(1,3,3);
Im_denoise = uint8(255*mat2gray(De_noise_im./Denominator));
imshow(Im_denoise);
xlabel(['PSNR = ', num2str(psnr(Im_denoise,original_im))])
title('Denoising Image');




end