function MPexample

close all
clear all
clc
%���� DCT dictionary
N = 1500; %DCT dictionary atoms ���Ӽ�
signal_length = 1024; %�T������
dictionary = zeros(signal_length, N);
dictionary(:, 1) = 1/(sqrt(N))*ones(signal_length,1);
for i = 2 : N
    dictionary(:, i) = sqrt(2/N)*cos((pi/N)*((1:signal_length)' + (1/2))*(i-1));
end

% ���J���հT��
load cuspamax;

% matching pursuit
[coe, a_atoms] = MP(cuspamax', dictionary, 32);

%ø��
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
    g = abs(D'*x);
    [val, ind] = max(g); %��X���n�ȳ̤j�����q
    n_val = val / (x(:,ind)'*x(:,ind)); %�p���v��
    a = [a, n_val]; %�s����v��
    a_atoms = [a_atoms, x(:,ind)];%�s����v�ȹ������@atom
    R = R - n_val*x(:,ind);
end
end