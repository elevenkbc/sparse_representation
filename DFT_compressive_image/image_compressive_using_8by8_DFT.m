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

mDFT8 = complex(zeros(8, 8)); %每一個column 為DCT 在 frequency space 的 basis

N8 = 8;
omega = exp((-2*pi*j)/N8);
for k = 1 : N8
    for i = 1 : N8 
        mDFT8(i, k) = (1/sqrt(N8))*omega^((i-1)*(k-1)); 
    end
end

imDFT8 = complex(zeros(8, 8)); %每一個column 為DCT 在 frequency space 的 basis
omega1 = exp((2*pi*j)/8);
for k = 1 : 8
    for i = 1 : 8 
        imDFT8(i, k) = (1/sqrt(8))*omega1^((i-1)*(k-1));
    end
end


%plot 8by8 2D DFT basis

DFTbasis_real = zeros(64,64);
DFTbasis_image = zeros(64,64);
for i = 1 : 8 %vertical index
    for j = 1 : 8 %horzontal index
        block_ind_start_i = (i-1)*8 + 1;
        block_ind_start_j = (j-1)*8 + 1;
        DFTbasis_real((block_ind_start_i:block_ind_start_i+7), (block_ind_start_j:block_ind_start_j+7)) = real(mDFT8(:,i)*mDFT8(:,j).');
        DFTbasis_image((block_ind_start_i:block_ind_start_i+7), (block_ind_start_j:block_ind_start_j+7)) = imag(mDFT8(:,i)*mDFT8(:,j).');
    end
end
plot_8by8_grid(DFTbasis_real);
title('FDT  2D basis real part');
plot_8by8_grid(DFTbasis_image);
title('FDT  2D basis image part');


DFT_coe = zeros(im_ver, im_hor);
compressive_im = zeros(im_ver, im_hor);

%對於每一個 8by8 Block 做2D DFT
complexA = complex(double(A)); %轉成 complex 才能計算 DFT
for i = 1 : block_ver_num
    for j = 1 : block_hor_num
        %提出8by8 block
        block_ind_start_i = (i-1)*8 + 1;
        block_ind_start_j = (j-1)*8 + 1;
        temp_block = complexA((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7));
        
        temp_dct = (mDFT8')*temp_block*mDFT8; %做 8by8 2D DFT
        DFT_coe((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7)) = temp_dct; %將係數塞回去
        % DCT 係數 threadhold
%         thrn = 7;
%         Thr_matrix = fliplr(toeplitz(zeros(1,8), [zeros(1, 8-thrn), ones(1,thrn)]));
%         
%         temp_dct_thr = temp_dct.*Thr_matrix;
        
        temp_dct_thr = temp_dct.*(abs(temp_dct) > 20);
        %IDCT
        
        compressive_im((block_ind_start_i:block_ind_start_i+7),(block_ind_start_j:block_ind_start_j+7)) = imDFT8'*temp_dct_thr*imDFT8;
    end
end
figure; imagesc(abs(DFT_coe)); colormap('gray')
title('DFT 係數')

figure; imagesc(abs(complexA)); colormap('gray')
title('原始圖片')

Ps = psnr(uint8(abs(compressive_im)), A);
figure; imagesc(abs(compressive_im)); colormap('gray')
title('compressive 圖片')
xlabel(['PSNR = ', num2str(Ps)]);


