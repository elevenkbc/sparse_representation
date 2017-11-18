
close all
clear all
clc

% 載入測試訊號
load cuspamax;

%產生 DCT dictionary
N = 1500; %DCT dictionary atoms 的個數
signal_length = length(cuspamax); %訊號長度
dictionary = zeros(signal_length, N);
for k = 1 : N
    if k == 1
        dictionary(:, k) = (sqrt(1/signal_length))*cos((pi/N)*(k-1)*((0:(signal_length-1))' + (1/2)*ones(signal_length,1)));
    else
        Temp = cos((pi/N)*(k-1)*((0:(signal_length-1))' + (1/2)*ones(signal_length,1)));
        Temp = Temp - mean(Temp);
        dictionary(:, k) = Temp/norm(Temp);
    end
end


% matching pursuit
L = 40; %OMP 中 0norm 的最大可能值
coe = OMP(dictionary, cuspamax',  L);

%繪圖
figure
plot(cuspamax,'b', 'linewidth', 1.5);
hold on
plot(dictionary*coe,'r--', 'linewidth', 1.5);

legend('original signal', 'OMP-by-DCT dictionary')

%計算 coe 中的非零元數
spark = 0;
for i = 1 : N
    if coe(i) ~= 0
        spark = spark +1;
    end
end
title(['||a||_0 = ', num2str(spark)]);
xlim([1, signal_length]);

