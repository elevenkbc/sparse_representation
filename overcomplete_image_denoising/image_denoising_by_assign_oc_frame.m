function image_denoising_by_assign_oc_frame(assign_frame, original_im, noising_im, Sigma)
d_noising_im = double(noising_im);
[im_ver, im_hor] = size(d_noising_im);
Block_im_n = im2col(d_noising_im, [8, 8],'sliding'); %�N�t�����T���v���A�Υi�H���|��8by8 block ���X�ӡA�Ʀ��x�}��column
%�C��Block_im_n �h��������
%�åB�N�C�� Block �������Ȧs�� Block_im_n_mean
Block_im_n_mean = mean(Block_im_n);
Block_im_n_mean0 = Block_im_n - kron(Block_im_n_mean, ones(64,1));


%% Denoising �{�Ƕ}�l!!
%���q1  ��overcomplete DCT �D�X�C�Ӥ����v���� sparse��ܫY��
Block_sparse_coe = OMPdenoise(assign_frame, Block_im_n_mean0, Sigma); %�o���n�]�ܤ[
 
%���q2  �b�w���C�ӧt���T�������x�}����ܫY�ƪ����p�U�A��close-form �D�Xdenoise �v�� hat_A
%�����M�� close-form ����
% close-form �ݰ_�ӫܽ����A�ƹ�W���u�O��Ҧ������|��block����̷�
lam = 0.3; %�Ѽ�

%Denominator �����|���v���C��pixel �n���W���v��
%�c�y�v��
Denominator = lam*ones(im_ver, im_hor);
for i = 1 : (im_ver-8+1)
    for j = 1 : (im_hor-8+1)
        Denominator(i:(i+7), j:(j+7)) = Denominator(i:(i+7), j:(j+7)) + ones(8,8);
    end
end
%�c�y�����|�� sparse reconstruction
Recover_block_IM = assign_frame*Block_sparse_coe + Block_im_n_mean;
De_noise_im = lam*d_noising_im;
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