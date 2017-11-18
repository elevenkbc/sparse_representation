

close all
clear all
clc

% 載入測試訊號
load cuspamax;


%產生 DCT dictionary
N = 1500; %DCT dictionary atoms 的個數
signal_length = length(cuspamax); %訊號長度
dictionary = zeros(signal_length, N);
dictionary(:, 1) = 1/(sqrt(N))*ones(signal_length,1);
for i = 2 : N
    dictionary(:, i) = sqrt(2/N)*cos((pi/N)*((0:(signal_length-1))' + (1/2))*(i-1));
end


% matching pursuit
L = 40;
[coe, a_atoms] = MP(cuspamax', dictionary, L);

%繪圖
figure
plot(cuspamax,'b', 'linewidth', 1.5);
hold on
plot(a_atoms*coe','r--', 'linewidth', 1.5);

legend('original signal', 'MP-by-DCT dictionary')
title(['||a||_0 = ', num2str(L)]);
xlim([1, signal_length]);

