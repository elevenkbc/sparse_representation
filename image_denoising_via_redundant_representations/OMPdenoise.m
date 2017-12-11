function [A] = OMPdenoise(D, X, Im_sigma)
% Orthogonal matching pursuit  (greed algorithm)   
%  
%   min_{a}  ||a||_0        subject to  || x_{ij} - D a ||_2 <= C*sigma^2
%   C�Px ���צ����Y
%
% input :
% X ��l�T����@column�Ʀ����T���x�} X = [x_1, x_2, x_3, ..., x_p];
% D �T���r�� (DCT, harr wavelet ����)
% L a���D�s���Ӽ�
% output :
% a ��ܫY��
% �p�G || x - D a ||_2 < C*sigma^2 �h���X�j��
waitBarOn = 1;
if (waitBarOn)
    counterForWaitBar = size(X,2);
    h = waitbar(0,'OMP In Process ...');
end
[n, p] = size(X);
C = n*1.2;
k = size(D, 2);
A = sparse(k, p);
L = ceil(n/2);
%���ND ���C�@��atom����
norm_D = zeros(size(D));
for i = 1 : k
    norm_D(:,i) = D(:,i)/norm(D(:,i));
end

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

        temp_coe = pinv(Phi)*X(:,i); %�p��N X(:,j) ��v�� R(Phi) �Ŷ��ϱo|| X(:,i) - Phi*(temp_coe)|| �̤p���Y��
        P = Phi*pinv(Phi); %P ����v�� R(Phi) �Ŷ��������v�x�}  P = Phi*inv(Phi'*Phi)*Phi'
        R = X(:,i) - P*X(:,i); %R ���ѤU���ݮt�V�q�AR�ݺ����P �w�g�ιL�� atom �����۫���
        
        if norm(R)^2 < C*Im_sigma^2 %�p�G�ѤU���ݮt�q�� norm ���� �p�� C sigma^2�A�h���X�j��
            break;
        end
    end
    %�N�Y�ƶ�^ �Y�Ưx�}��
    A(atom_ind, i) = temp_coe;
    if (waitBarOn)
        waitbar(i/counterForWaitBar);
    end
end
if (waitBarOn)
    close(h);
end
end