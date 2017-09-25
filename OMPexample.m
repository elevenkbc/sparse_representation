function OMPexample

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
[coe, a_atoms] = OMP(cuspamax', dictionary, 35);

%ø��
figure
plot(cuspamax,'b', 'linewidth', 1.5);
hold on
plot(a_atoms*coe','r--', 'linewidth', 1.5);

legend('original signal', 'OMP-by-DCT dictionary')
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
    g = abs(D'*x);
    [val, ind] = max(g); %��X���n�ȳ̤j�����q
    n_val = val / (x(:,ind)'*x(:,ind)); %�p���v��
    a_atoms = [a_atoms, x(:,ind)];%�s����v�ȹ������@atom
    a = [a, n_val]; %�s����v��
    Phi = [Phi, D(:,ind)]; %Phi �Ψ��x�s�A���e�ҿ�ܹL�� atom

    P = Phi*pinv(Phi); %P ����v�� Phi ��Ŷ��������v�x�}  P = Phi*inv(Phi'*Phi)*Phi'
    
    R = (eye(x_len) - P)*R;
end
end