function [Signal] = iwavelet_tf_Morlet(TF_coe, signal_t, layer)
Signal = zeros(1, length(signal_t));
for j = 1 : layer
    for n = 0 : (2^(j-1)-1)
        Morlet_psi_j_n = (sqrt(2^(j-1)))*Morlet(2^(j-1)*signal_t - n);
        Signal = Signal + TF_coe(1, 2^(j-1)+n)*Morlet_psi_j_n;
    end
end
end
