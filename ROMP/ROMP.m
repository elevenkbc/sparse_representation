function [A] = ROMP(D, X, L)
% Regularized orthogonal matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% X 原始訊號當作column排成的訊號矩陣 X = [x_1, x_2, x_3, ..., x_p];
% D 訊號字典 (DCT, harr wavelet 等等)
% L a的非零元個數
% output :
% a 表示係數

[n, p] = size(X);
k = size(D, 2);
A = zeros(k, p);

%先將D 中每一個atom單位化
norm_D = zeros(size(D));
for i = 1 : k
    norm_D(:,i) = D(:,i)/norm(D(:,i));
end

for i = 1 : p
    %initailize
    R = X(:, i);
    atom_ind = [];%初始下標集合為空集合
    for j = 1 : L
        g = abs(norm_D'*R); %計算內積絕對值
        if (sum(g~=0) <= L) %如果g的非零元小於L，則選取全部的下標
            t1 = (1:size(g,1)).*(g~=0); 
            current_index = t1(t1~=0);
        else
            %找出內積絕對值最大的L個分量，並且將L個分量對應的下標存成 max_L_index
            [sorted_g , total_index]  =sort(g ,'descend');
            max_L_index = total_index(1:L);
            current_index = max_L_index;
        end
        
        R_current_index = Regular_index(current_index, g); %將下標正則化，這一步會搜尋所有的power set，故L不可以太大
        
        if length([atom_ind, R_current_index'])>= 2*L %如果 |atom_ind|>= 2L，跳出迴圈
            break;
        else
            atom_ind = [atom_ind, R_current_index']; %更新下標集
            D_atom_ind = norm_D(:,atom_ind); %更新下標集所對應的字典
            
            PinvD = pinv(D_atom_ind);
            temp_coe = PinvD*X(:,i); %計算將 X(:,j) 投影到 R(D_atom_ind) 空間使得|| X(:,i) - Phi*(temp_coe)|| 最小的係數
            P = D_atom_ind*PinvD; %P 為投影到 R(D_atom_ind) 空間的正交投影矩陣，投影矩陣的公式為 P = D*inv(D'*D)*D'
            R = X(:,i) - P*X(:,i); %R 為剩下的殘差向量，R需滿足與 D_atom_ind 每個column均互相垂直
        end
    end
    %將係數填回 係數矩陣Ａ
    A(atom_ind, i) = temp_coe;
end
end
function [J0]= Regular_index(J, U)
%J0 為正則化後的J座標子集
% 正則化的座標子集 J0 為以下最佳化問題的解
% J0 = arg max { sum_{j in J0} |U_j|^2 }
% subject to  |u(i)| <= 2 |u(j)| for all i, j in J0,  where J0 in P(J)
% P(J) 表示J的power set，為集合J的所有可能子集的收集起來的集合

%其中 J, U 都是行向量
length_J = length(J);
B = ff2n(length_J); %B的每一個row 為 J的一種子集可能
max_norm_value = 0;
test_value_array = zeros(1, size(B,1));
for i = 1 : size(B, 1)
    temp_index = J(B(i,:) == 1);
    sub_U = U(temp_index);
    if ~isempty(sub_U)
        if (2*min(abs(sub_U)) >= max(abs(sub_U)))&&(sum(sub_U.*sub_U) > max_norm_value)
            max_norm_value = sum(sub_U.*sub_U);
            J0 = temp_index;
            test_value_array(i) = sum(sub_U.*sub_U);
        end
    end
end
end