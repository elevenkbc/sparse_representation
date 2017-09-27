function OMP_image_sparse_representation
close all
clear all
clc

%Sparce representation donoising
% imput noising data
A = imread('horizontal_test.bmp');
double_A = double(A); %�Nuint8 �ର Double
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


%�N2D DCT �� atoms �H��V�q�Ʀ��x�}
DCT_Dictionary = zeros(M*N, P*Q);

for i = 1 : P*Q
    Temp = DCT_2D(:, :, i);
    DCT_Dictionary(:,i) = Temp(:);   
end


%�}�l MP denoising by 2D-DCT frame
A_col = double_A(:); %�N�v���Ʀ���V�q

L = 50;
[coe, a_atoms] = OMP(A_col, DCT_Dictionary, L);
denoising_im_col = a_atoms*coe';
denoising_im = reshape(denoising_im_col, M, N);

%ø��
figure;
subplot(1,2,1)
imshow(A); 
title('��l�Ϥ�')
subplot(1,2,2)
imshow(uint8(denoising_im));
P1 = psnr(A, uint8(denoising_im));
title(['OMP by DCT with ||a||_0 = ', num2str(L)]);
xlabel(['PSNR = ', num2str(P1)])

end

function output = scale_factor(input, M)
if input == 0
    output = 1/sqrt(M);
else
    output = sqrt(2/M);
end
end

function [a, a_atoms] = OMP(x, D, L)
% Orthogonal matching pursuit  (greed algorithm)   
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
Phi = [];

x_len = length(x);
for i = 1 : L
    g = abs(D'*R);
    [val, ind] = max(g); %��X���n�ȳ̤j�����q
    n_val = (R'*D(:,ind)) / (D(:,ind)'*D(:,ind)); %�p���v��
    a = [a, n_val]; %�s����v��
    a_atoms = [a_atoms, D(:,ind)];%�s����v�ȹ������@atom
    Phi = [Phi, D(:,ind)]; %Phi �Ψ��x�s�A���e�ҿ�ܹL�� atom

    P = Phi*pinv(Phi); %P ����v�� Phi ��Ŷ��������v�x�}  P = Phi*inv(Phi'*Phi)*Phi'
    
    R = (eye(x_len) - P)*R;
end
end