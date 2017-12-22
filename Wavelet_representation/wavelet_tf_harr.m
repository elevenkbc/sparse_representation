function [TF_coe] = wavelet_tf_harr(signal_t, signal, layer)
%�����Y�� �s��wavelet �ഫ
TF_coe = [];
for j = 1 : layer
    for n = 0 : (2^(j-1)-1)
        Harr_psi_j_n = sqrt(2^(j-1))*Harr(2^(j-1)*signal_t - n);
        TF_coe = [TF_coe, trapz(signal_t, signal.*Harr_psi_j_n)];
    end
end

end
