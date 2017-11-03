function [a, a_atoms] = MP(x, D, L)
% Matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% x 原始訊號
% D 訊號字典 (DCT, harr wavelet 等等)
% L a的非零元個數
% output :
% a 表示係數
% a_atoms 與係數相關的字典中的"詞"


a = [];
a_atoms = [];

%initailize

R = x;

for i = 1 : L
    g = abs(D'*R);
    [val, ind] = max(g); %找出內積值最大的分量
    n_val = (R'*D(:,ind)) / (D(:,ind)'*D(:,ind)); %計算投影值
    a = [a, n_val]; %存取投影值
    a_atoms = [a_atoms, D(:,ind)];%存取投影值對應的　atom
    R = R - n_val*D(:,ind);
end
end