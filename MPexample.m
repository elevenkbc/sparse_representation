function MPexample

close all
clear all
clc
%產生 DCT dictionary
N = 1500; %DCT dictionary atoms 的個數
signal_length = 1024; %訊號長度
dictionary = zeros(signal_length, N);
dictionary(:, 1) = 1/(sqrt(N))*ones(signal_length,1);
for i = 2 : N
    dictionary(:, i) = sqrt(2/N)*cos((pi/N)*((1:signal_length)' + (1/2))*(i-1));
end

% 載入測試訊號
load cuspamax;

% matching pursuit
[coe, a_atoms] = MP(cuspamax', dictionary, 32);

%繪圖
figure
plot(cuspamax,'b', 'linewidth', 1.5);
hold on
plot(a_atoms*coe','r--', 'linewidth', 1.5);

legend('original signal', 'MP-by-DCT dictionary')
end
function [a, a_atoms] = MP(x, D, L)
% Matching pursuit  (greed algorithm)   
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

for i = 1 : L
    g = abs(D'*x);
    [val, ind] = max(g); %找出內積值最大的分量
    n_val = val / (x(:,ind)'*x(:,ind)); %計算投影值
    a = [a, n_val]; %存取投影值
    a_atoms = [a_atoms, x(:,ind)];%存取投影值對應的　atom
    R = R - n_val*x(:,ind);
end
end