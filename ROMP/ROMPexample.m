
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


% Regularized orthogonal matching pursuit
L = 6; %���T�� sparsity �����p�ȡA�`�N���i�H�Ӥj�A�]���bRegularized index ���L�{���A�|�j�M�Ҧ�power set
coe = ROMP(dictionary, cuspamax',  L);

%ø��
figure
plot(cuspamax,'b', 'linewidth', 1.5);
hold on
plot(dictionary*coe,'r--', 'linewidth', 1.5);

legend('Original signal', 'Recovered signal from ROMP', 'location', 'Best')
%�p�� coe �����D�s����
spark = 0;
for i = 1 : N
    if coe(i) ~= 0
        spark = spark +1;
    end
end
title(['||a||_0 = ', num2str(spark)]);
xlabel(['MSE = ', num2str(sum((cuspamax'-dictionary*coe).^2)/length(cuspamax))]);
xlim([1, signal_length]);



