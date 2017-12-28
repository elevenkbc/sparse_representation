function [A] = OMP(D, X, L)
% Orthogonal matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% X ��l�T����@column�Ʀ����T���x�} X = [x_1, x_2, x_3, ..., x_p];
% D �T���r�� (DCT, harr wavelet ����)�A�r�媺�C�@��column��norm�����O1�C
% L a���D�s���Ӽ�
% output :
% a ��ܫY��
% �p�G || x - D a ||_2 < tol �h���X�j��
[n, p] = size(X);
k = size(D, 2);
A = zeros(k, p);


for i = 1 : p
    %initailize
    R = X(:, i);
    Phi = [];
    atom_ind = [];
    for j = 1 : L
        g = abs(D'*R);
        [val, ind] = max(g); %��X���n����ȳ̤j�����q
        atom_ind = [atom_ind, ind];% atom_ind ��ܦb��i���|�N�H�e�A�Ҧ��w�g�ιL�� atom index�C
        Phi = [Phi, D(:,ind)]; %Phi �Ψ��x�s�A���e�ҿ�ܹL�� atom
        
        %%% �N X(:, j) ���Ѧ��P R(Phi) �Ŷ� �P ���� R(Phi)���Ŷ� �����M(direct sum)�C(Phi �x�}��Range space �O�� R(Phi))
        PinvPhi = pinv(Phi);
        temp_coe = PinvPhi*X(:,i); %�p��N X(:,j) ��v�� R(Phi) �Ŷ��ϱo|| X(:,i) - Phi*(temp_coe)|| �̤p���Y��
        P = Phi*PinvPhi; %P ����v�� R(Phi) �Ŷ��������v�x�}  P = Phi*inv(Phi'*Phi)*Phi'
        R = X(:,i) - P*X(:,i); %R ���ѤU���ݮt�V�q�AR�ݺ����P �w�g�ιL�� atom �����۫���
        
        if norm(R) < 1.0e-6 %�p�G�ѤU���ݮt�q�Ӥp�A�N���X�j��
            break;
        end
    end
    %�N�Y�ƶ�^ �Y�Ưx�}��
    A(atom_ind, i) = temp_coe;
end
end