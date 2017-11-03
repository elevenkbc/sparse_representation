function OMPexample
close all
clear all
clc

% ���J���հT��
load cuspamax;

%���� DCT dictionary
N = 1500; %DCT dictionary atoms ���Ӽ�
signal_length = length(cuspamax); %�T������
dictionary = zeros(signal_length, N);
dictionary(:, 1) = 1/(sqrt(N))*ones(signal_length,1);
for i = 2 : N
    dictionary(:, i) = sqrt(2/N)*cos((pi/N)*((0:(signal_length-1))' + (1/2))*(i-1));
end


% matching pursuit
L = 40; %OMP �� 0norm ���̤j�i���
[coe, a_atoms] = OMP(cuspamax', dictionary, L, 1.0e-1);

%ø��
figure
plot(cuspamax,'b', 'linewidth', 1.5);
hold on
plot(a_atoms*coe','r--', 'linewidth', 1.5);

legend('original signal', 'OMP-by-DCT dictionary')
title(['||a||_0 = ', num2str(length(coe))]);
xlim([1, signal_length]);
end
