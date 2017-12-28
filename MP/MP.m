function [A] = MP(D, X, L)
% Matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% x ��l�T��
% D �T���r�� (DCT, harr wavelet ����)�AD���C�@��column��norm�����O1
% L a���D�s���Ӽ�
% output :
% a ��ܫY��
% a_atoms �P�Y�Ƭ������r�夤��"��"


A = sparse(zeros(size(D,2), size(X,2)));


%initailize
for j = 1 : size(X, 2)
    R = X(:,j);
    for i = 1 : L
        g = abs((D')*R);
        [val, ind] = max(g); %��X���n�ȳ̤j�����q
        A(ind, j) = A(ind, j) + (D(:,ind)')*R; %�s����v��
        R = R - A(ind, j)*D(:,ind);
        if norm(R) < 1.0e-6 %�p�G�ѤU���ݮt�q�Ӥp�A�N���X�j��
            break;
        end
    end
end
end