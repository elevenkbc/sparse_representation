function [A] = MP(D, X, L)
% Matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% x 原始訊號
% D 訊號字典 (DCT, harr wavelet 等等)，D的每一個column的norm必須是1
% L a的非零元個數
% output :
% a 表示係數
% a_atoms 與係數相關的字典中的"詞"


A = sparse(zeros(size(D,2), size(X,2)));


%initailize
for j = 1 : size(X, 2)
    R = X(:,j);
    for i = 1 : L
        g = abs((D')*R);
        [val, ind] = max(g); %找出內積值最大的分量
        A(ind, j) = A(ind, j) + (D(:,ind)')*R; %存取投影值
        R = R - A(ind, j)*D(:,ind);
        if norm(R) < 1.0e-6 %如果剩下的殘差量太小，就跳出迴圈
            break;
        end
    end
end
end