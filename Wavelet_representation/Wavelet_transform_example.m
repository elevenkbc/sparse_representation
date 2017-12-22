close all
clear all
clc
%產生test signal
t = 0:0.001:1;
x = sin(2*pi*2*t);
t1 = 300;
t2 = 700;
x(t1:t2) = x(t1:t2) - 0.3;
x = x + 0.02*randn(size(x));

layer = 8; %Wavelet 最大的層數 在地5層有 2^5個係數

[Coe_harr] = wavelet_tf_harr(t, x, layer); %對x做CWT，以morlet為小波基底
[TF_harr] = Wavelet_coe_TF(abs(Coe_harr), layer); %輸出TF 表示圖
[Recovered_signal_harr] = iwavelet_tf_harr(Coe_harr, t, layer); %由Wavelet係數，與harr基底近似原本訊號

[Coe_morlet] = wavelet_tf_Morlet(t, x, layer);  %對x做CWT，以morlet為小波基底
[TF_morlet] = Wavelet_coe_TF(abs(Coe_morlet), layer);%輸出TF 表示圖
[Recovered_signal_morlet] = iwavelet_tf_Morlet(Coe_morlet, t, layer); %由Wavelet係數，與Morlet基底近似原本訊號


%開始繪圖
figure('pos',[220 80 1000 950]); 
subplot(3,2,1); plot(t, real(x)); title('original signal'); axis tight;
subplot(3,2,2); plot(t, real(x)); title('original signal'); axis tight;
subplot(3,2,3); imagesc(TF_harr); ylabel('scale'); title('TF form Harr wavelet'); subplot(3,2,4); imagesc(TF_morlet); ylabel('scale');  title('TF form Morlet wavelet');
subplot(3,2,5); plot(t, Recovered_signal_harr); title('reconstruction signal by Harr wavelet');
subplot(3,2,6); plot(t, real(Recovered_signal_morlet), t, imag(Recovered_signal_morlet)); title('reconstruction signal by Morlet wavelet'); legend('real part', 'imagine part', 'Location', 'Best');

