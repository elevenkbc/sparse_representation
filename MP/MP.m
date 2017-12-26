function [A] = MP(D, X, L)
% Matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% x ﹍癟腹
% D 癟腹ㄥ (DCT, harr wavelet 单单)D–columnnormゲ斗琌1
% L a獶箂じ计
% output :
% a ボ玒计
% a_atoms 籔玒计闽ㄥい"迭"


A = sparse(zeros(size(D,2), size(X,2)));


%initailize
for j = 1 : size(X, 2)
    R = X(:,j);
    for i = 1 : L
        g = abs((D')*R);
        [val, ind] = max(g); %тず縩程だ秖
        A(ind, j) = A(ind, j) + (D(:,ind)')*R; %щ紇
        R = R - A(ind, j)*D(:,ind);
    end
end
end