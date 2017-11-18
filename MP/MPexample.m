

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
L = 40;
[coe, a_atoms] = MP(cuspamax', dictionary, L);

%ø��
figure
plot(cuspamax,'b', 'linewidth', 1.5);
hold on
plot(a_atoms*coe','r--', 'linewidth', 1.5);

legend('original signal', 'MP-by-DCT dictionary')
title(['||a||_0 = ', num2str(L)]);
xlim([1, signal_length]);

