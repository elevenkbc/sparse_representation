function MP_image_sparse_denoising

close all
clear all
clc

%Sparce representation donoising
% imput noising data
A = imread('lena_gray_86.bmp');
double_A = double(A); %�Nuint8 �ର Double
double_A = double(uint8(double_A + 20*randn(size(double_A))));
[M, N] = size(double_A);


%����2D DCT-II dictionary
P = 150; Q = 150; %�`�@���� P*Q �� DCT-II Frame
DCT_2D = zeros(M, N, P*Q);
tic
for p = 1 : P
    for q = 1 : Q    
        DCT_2D(:, :, p + (q-1)*Q) = scale_factor((p-1),P)*scale_factor((q-1),Q)*...
            cos((pi*(2*(0:(M-1))'+1)*(p-1))./(2*M))*cos((pi*(2*(0:(N-1))+1)*(q-1))/(2*N));
    end
end
toc

ZZThr = 1000;
ZZmatrix = reshape(1:(P*Q),P,Q);
ZZindex = zigzag(ZZmatrix);
ZZindex_cut = ZZindex(1:ZZThr);



%�N2D DCT �� atoms �H��V�q�Ʀ��x�}
DCT_Dictionary = zeros(M*N, P*Q);

for i = 1 : P*Q
    Temp = DCT_2D(:, :, i);
    DCT_Dictionary(:,i) = Temp(:);   
end
DCT_Dictionary_cut = DCT_Dictionary(:, ZZindex_cut);


%�}�l MP denoising by 2D-DCT frame
A_col = double_A(:); %�N�v���Ʀ���V�q

L = 1000;
[coe, a_atoms] = MP(A_col, DCT_Dictionary, L);
denoising_im_col = a_atoms*coe';
denoising_im = reshape(denoising_im_col, M, N);

%ø��
figure;
subplot(1,3,1)
imshow(uint8(A)); 
title('��l�Ϥ�');
subplot(1,3,2)
imshow(uint8(double_A)); 
P2 = psnr(uint8(double_A), A);
title('noising �Ϥ�');
xlabel(['PSNR = ', num2str(P2)]);

subplot(1,3,3)
imshow(uint8(denoising_im));
P2 = psnr(uint8(denoising_im), A);
title(['MP by DCT with ||a||_0 = ', num2str(L)]);
xlabel(['PSNR = ', num2str(P2)]);

end

function output = scale_factor(input, M)
if input == 0
    output = 1/sqrt(M);
else
    output = sqrt(2/M);
end
end

function [a, a_atoms] = MP(x, D, L)
% Matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% x ��l�T��
% D �T���r�� (DCT, harr wavelet ����)
% L a���D�s���Ӽ�
% output :
% a ��ܫY��
% a_atoms �P�Y�Ƭ������r�夤��"��"


a = [];
a_atoms = [];

%initailize

R = x;

for i = 1 : L
    g = abs(D'*R);
    [val, ind] = max(g); %��X���n�ȳ̤j�����q
    n_val = (R'*D(:,ind)) / (D(:,ind)'*D(:,ind)); %�p���v��
    a = [a, n_val]; %�s����v��
    a_atoms = [a_atoms, D(:,ind)];%�s����v�ȹ������@atom
    R = R - n_val*D(:,ind);
end
end