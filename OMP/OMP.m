function [a, a_atoms] = OMP(x, D, L, tol)
% Orthogonal matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% x ��l�T��
% D �T���r�� (DCT, harr wavelet ����)
% L a���D�s���Ӽ�
% output :
% a ��ܫY��
% a_atoms �P�Y�Ƭ������r�夤��"��"
% �p�G || x - D a ||_2 < tol �h���X�j��


a = [];
a_atoms = [];

%initailize

R = x;
Phi = [];

x_len = length(x); 
for i = 1 : L
    g = abs(D'*R);
    [val, ind] = max(g); %��X���n�ȳ̤j�����q
    n_val = (R'*D(:,ind)) / (D(:,ind)'*D(:,ind)); %�p���v��
    a = [a, n_val]; %�s����v��
    a_atoms = [a_atoms, D(:,ind)];%�s����v�ȹ������@atom
    if norm(x - a_atoms*(a')) < tol
        break;
    end
    Phi = [Phi, D(:,ind)]; %Phi �Ψ��x�s�A���e�ҿ�ܹL�� atom
    P = Phi*pinv(Phi); %P ����v�� Phi ��Ŷ��������v�x�}  P = Phi*inv(Phi'*Phi)*Phi'
    
    
    R = (eye(x_len) - P)*R;
end
end