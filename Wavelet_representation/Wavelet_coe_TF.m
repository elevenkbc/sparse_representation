function [TF_im] = Wavelet_coe_TF(TF_coe, layer)
TF_im = zeros(layer, 2^(layer));
for j = 1 : layer
    Temp_coe = TF_coe(2^(j-1):((2^j)-1));
    TF_im(j, :) = kron(Temp_coe, ones(1, 2^(layer - j + 1)));
end
end