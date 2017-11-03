clear all;
close all;
clc
A = imread('barbara.bmp');
A_n = imread('barbara_n.bmp');
[im_ver, im_hor] = size(A);
block_ver_num = im_ver/8;
block_hor_num = im_hor/8;



% 以 DCT-II overcomplete dictionary
% 利用 OMP 做 sparse coding 來找出係數與相關的 atom
% sparse representation 可以用來做　image denoising
% 注意 DCT-II overcomplete dictionary 在　barbara denoising 效果比較好


% 構造overcomplete DCT frame
N_1 = 8;
N_2 = 16;
% DCT-II
mDCT_oc = zeros(N_1, N_2);
for k = 1 : N_2
    if k == 1
        mDCT_oc(:, k) = (sqrt(1/N_1))*cos((pi/N_2)*(k-1)*((0:(N_1-1))' + (1/2)*ones(N_1,1)));
    else
        mDCT_oc(:, k) = (sqrt(2/N_1))*cos((pi/N_2)*(k-1)*((0:(N_1-1))' + (1/2)*ones(N_1,1)));
    end
    
end


%plot 8by8 2D DCT basis
DCT_2D_frame = zeros(N_1*N_2,N_1*N_2);
for i = 1 : N_2 %vertical index
    for j = 1 : N_2 %horzontal index
        block_ind_start_i = (i-1)*8 + 1;
        block_ind_start_j = (j-1)*8 + 1;
        DCT_2D_frame((block_ind_start_i:block_ind_start_i+7), (block_ind_start_j:block_ind_start_j+7)) = mDCT_oc(:,i)*mDCT_oc(:,j)';
    end
end
plot_8by8_grid(DCT_2D_frame)


%對於每一個 8by8 Block 做 sparse coding

%create codebook D
D = zeros(N_1*N_1, N_2);
for i = 1 : N_2 %vertical index
    for j = 1 : N_2 %horzontal index
        temp_2D = mDCT_oc(:,i)*mDCT_oc(:,j)';
        D(:, (i-1)*N_2 + j) = temp_2D(:);
    end
end
doubleA = double(A); %轉成 double 才能計算 OMP


noise_imag = double(A_n);
denoising_im = zeros(im_ver, im_hor);

for i = 1 : block_ver_num
    for j = 1 : block_hor_num
        %提出8by8 block
        block_ind_start_i = (i-1)*8 + 1;
        block_ind_start_j = (j-1)*8 + 1;
        temp_block = noise_imag((block_ind_start_i:block_ind_start_i+(N_1 - 1)),(block_ind_start_j:block_ind_start_j+(N_1 - 1)));
        ob_col = temp_block(:);
        
        %sparse coding  using OMP
        L = N_2; %||coe||_0  <= L
        [coe, corres_D] = OMP(ob_col, D, L, 90.89); % 90.89 是denoising 效果比較好的參數
                
        
        %sparse representation
        S_R = corres_D*(coe');
        %將S_R 排回去 8 by 8
        S_R_8by8 = reshape(S_R, 8, 8);
        denoising_im((block_ind_start_i:block_ind_start_i+(N_1 - 1)),(block_ind_start_j:block_ind_start_j+(N_1 - 1))) = S_R_8by8;
        
        
    end
end




figure('Position',[100,300,1100,350]) 
subplot(1,3,1); imshow(A); title('原始圖片');
subplot(1,3,2); imshow(A_n); title('含有雜訊圖片');
Ps1 = psnr(A_n, A);
xlabel(['PSNR = ', num2str(Ps1)]);
subplot(1,3,3); imshow(uint8(denoising_im)); title('denoising by overcomplete DCT frame');
Ps2 = psnr(uint8(denoising_im), A);
xlabel(['PSNR = ', num2str(Ps2)]);




