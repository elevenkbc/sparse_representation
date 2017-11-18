
close all
clear all
clc

% ���J���հT��
load cuspamax;

%���� DCT dictionary
N = 1500; %DCT dictionary atoms ���Ӽ�
signal_length = length(cuspamax); %�T������
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
L = 40; %OMP �� 0norm ���̤j�i���
coe = OMP(dictionary, cuspamax',  L);

%ø��
figure
plot(cuspamax,'b', 'linewidth', 1.5);
hold on
plot(dictionary*coe,'r--', 'linewidth', 1.5);

legend('original signal', 'OMP-by-DCT dictionary')

%�p�� coe �����D�s����
spark = 0;
for i = 1 : N
    if coe(i) ~= 0
        spark = spark +1;
    end
end
title(['||a||_0 = ', num2str(spark)]);
xlim([1, signal_length]);

