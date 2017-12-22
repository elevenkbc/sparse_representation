function [TF_coe] = wavelet_tf_Morlet(signal_t, signal, layer)
%離散係數 連續wavelet 轉換  using Morlet
TF_coe = [];
for j = 1 : layer
    for n = 0 : (2^(j-1)-1)
        Morlet_psi_j_n = (sqrt(2^(j-1)))*Morlet(2^(j-1)*signal_t - n);
        TF_coe = [TF_coe, trapz(signal_t, signal.*conj(Morlet_psi_j_n))];
    end
end

end