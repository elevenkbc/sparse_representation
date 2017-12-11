clear all;
close all;
clc
A = imread('boat.bmp');
[im_ver, im_hor] = size(A);

% �H DCT-II overcomplete dictionary
% �Q�� OMP �� sparse coding �ӧ�X�Y�ƻP������ atom

% �c�yovercomplete DCT frame
N_1 = 8;  %�C��Block ��8by8
N_2 = 16;  %Dictionary �� 16*16 �� 8by8 2D-DCT �զ� 
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

% �s�y DCT_2D overcomplete DCT dictionray
DCT_2D_frame = zeros(N_1*N_1, N_2*N_2);
for i = 1 : N_2 %vertical index
    for j = 1 : N_2 %horzontal index
        if j == 1
            Temp2D = mDCT_oc(:,i)*mDCT_oc(:,j)';
            DCT_2D_frame(:, (i-1)*N_2 + j) = Temp2D(:);
        else
        Temp2D = mDCT_oc(:,i)*mDCT_oc(:,j)';
        Temp_col = Temp2D(:) - mean(Temp2D(:))*ones(size(Temp2D(:))); %���C��8by8�������Ȭ�0
        DCT_2D_frame(:, (i-1)*N_2 + j) = Temp_col/norm(Temp_col); %normalizing
        end
    end
end
% DCT_2D overcomplete DCT dictionary �c�y����

Sigma = 20;
A_n = double(A) + Sigma*randn(im_ver, im_hor); %���ͧt�����T���v��
Block_im_n = im2col(A_n, [8 8],'sliding'); %�N�t�����T���v���A�Υi�H���|��8by8 block ���X�ӡA�Ʀ��x�}��column

%�C��Block_im_n demean
Block_im_n_mean = mean(Block_im_n);
Block_im_n_mean0 = Block_im_n - kron(Block_im_n_mean, ones(64,1));


%% Denoising �{�Ƕ}�l!!
%���q1  ��overcomplete DCT �D�X�C�Ӥ����v���� sparse��ܫY��
Block_sparse_coe = OMPdenoise(DCT_2D_frame, Block_im_n_mean0, Sigma); %�o���n�]�ܤ[
 
%���q2  �b�w���C�ӧt���T�������x�}����ܫY�ƪ����p�U�A��close-form �D�Xdenoise �v�� hat_A
%�����M�� close-form ����
% close-form �ݰ_�ӫܽ����A�ƹ�W���u�O��Ҧ������|��block����̷�
lam = 0.1; %�Ѽ�

%Denominator �����|���v���C��pixel �n���W���v��
%�c�y�v��
Denominator = lam*ones(im_ver, im_hor);
for i = 1 : (im_ver-8+1)
    for j = 1 : (im_hor-8+1)
        Denominator(i:(i+7), j:(j+7)) = Denominator(i:(i+7), j:(j+7)) + ones(8,8);
    end
end
%�c�y�����|�� sparse reconstruction
Recover_block_IM = DCT_2D_frame*Block_sparse_coe + Block_im_n_mean;
De_noise_im = lam*A_n;
coe_ind = 1; % Block_sparse_coe �� index
for i = 1 : (im_ver-8+1)
    for j = 1 : (im_hor-8+1)
        De_noise_im(j:(j+7), i:(i+7)) = De_noise_im(j:(j+7), i:(i+7)) + reshape(Recover_block_IM(:, coe_ind), [8,8]);
        coe_ind = coe_ind + 1; 
    end
end

%% �q�X���G
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


