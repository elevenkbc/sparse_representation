function [a, a_atoms] = MP(x, D, L)
% Matching pursuit  (greed algorithm)   
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


a = [];
a_atoms = [];

%initailize

R = x;

for i = 1 : L
    g = abs(D'*R);
    [val, ind] = max(g); %��X���n�ȳ̤j�����q
    n_val = (R'*D(:,ind)) / (D(:,ind)'*D(:,ind)); %�p���v��
    a = [a, n_val]; %�s����v��
    a_atoms = [a_atoms, D(:,ind)];%�s����v�ȹ������@atom
    R = R - n_val*D(:,ind);
end
end