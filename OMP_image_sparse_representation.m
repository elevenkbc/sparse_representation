function OMP_image_sparse_representation
close all
clear all
clc

%Sparce representation donoising
% imput noising data
A = imread('horizontal_test.bmp');
double_A = double(A); %將uint8 轉為 Double
[M, N] = size(double_A);


%產生2D DCT-II dictionary
P = 150; Q = 150; %總共產生 P*Q 個 DCT-II Frame
DCT_2D = zeros(M, N, P*Q);
tic
for p = 1 : P
    for q = 1 : Q    
        DCT_2D(:, :, p + (q-1)*Q) = scale_factor((p-1),P)*scale_factor((q-1),Q)*...
            cos((pi*(2*(0:(M-1))'+1)*(p-1))./(2*M))*cos((pi*(2*(0:(N-1))+1)*(q-1))/(2*N));
    end
end
toc


%將2D DCT 的 atoms 以行向量排成矩陣
DCT_Dictionary = zeros(M*N, P*Q);

for i = 1 : P*Q
    Temp = DCT_2D(:, :, i);
    DCT_Dictionary(:,i) = Temp(:);   
end


%開始 MP denoising by 2D-DCT frame
A_col = double_A(:); %將影像排成行向量

L = 50;
[coe, a_atoms] = OMP(A_col, DCT_Dictionary, L);
denoising_im_col = a_atoms*coe';
denoising_im = reshape(denoising_im_col, M, N);

%繪圖
figure;
subplot(1,2,1)
imshow(A); 
title('原始圖片')
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
% x 原始訊號
% D 訊號字典 (DCT, harr wavelet 等等)
% L a的非零元個數
% output :
% a 表示係數
% a_atoms 與係數相關的字典中的"詞"


a = [];
a_atoms = [];

%initailize

R = x;
Phi = [];

x_len = length(x);
for i = 1 : L
    g = abs(D'*R);
    [val, ind] = max(g); %找出內積值最大的分量
    n_val = (R'*D(:,ind)) / (D(:,ind)'*D(:,ind)); %計算投影值
    a = [a, n_val]; %存取投影值
    a_atoms = [a_atoms, D(:,ind)];%存取投影值對應的　atom
    Phi = [Phi, D(:,ind)]; %Phi 用來儲存，先前所選擇過的 atom

    P = Phi*pinv(Phi); %P 為投影到 Phi 行空間的正交投影矩陣  P = Phi*inv(Phi'*Phi)*Phi'
    
    R = (eye(x_len) - P)*R;
end
end