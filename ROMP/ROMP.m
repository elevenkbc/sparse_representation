function [A] = ROMP(D, X, L)
% Regularized orthogonal matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% X 原始訊號當作column排成的訊號矩陣 X = [x_1, x_2, x_3, ..., x_p];
% D 訊號字典 (DCT, harr wavelet 等等) 每個column必須長度為1
% L a的非零元個數
% output :
% a 表示係數

[n, p] = size(X);
k = size(D, 2);
A = zeros(k, p);


for i = 1 : p
    %initailize
    R = X(:, i);
    atom_ind = [];%初始下標集合為空集合
    for j = 1 : L    
        R_current_index = Regular_index(D'*R, L); %計算正則化的座標
        
        if length([atom_ind, R_current_index'])>= 2*L %如果 |atom_ind|>= 2L，跳出迴圈
            break;
        else
            atom_ind = [atom_ind, R_current_index']; %更新下標集
            D_atom_ind = D(:,atom_ind); %更新下標集所對應的字典
            
            PinvD = pinv(D_atom_ind);
            temp_coe = PinvD*X(:,i); %計算將 X(:,j) 投影到 R(D_atom_ind) 空間使得|| X(:,i) - Phi*(temp_coe)|| 最小的係數
            P = D_atom_ind*PinvD; %P 為投影到 R(D_atom_ind) 空間的正交投影矩陣，投影矩陣的公式為 P = D*inv(D'*D)*D'
            R = X(:,i) - P*X(:,i); %R 為剩下的殘差向量，R需滿足與 D_atom_ind 每個column均互相垂直
        end
        if norm(R) < 1.0e-6 %如果剩下的殘差量太小，就跳出迴圈
            break;
        end
    end
    %將係數填回 係數矩陣Ａ
    A(atom_ind, i) = temp_coe;
end
end
function [J0]= Regular_index(U, Kin)
%J0 為正則化後的J座標子集
% 正則化的座標子集 J0 為以下最佳化問題的解
% J0 = arg max { sum_{j in J0} |U_j|^2 }
% subject to  |u(i)| <= 2 |u(j)| for all i, j in J0,  where J0 in P(J)
% P(J) 表示J的power set，為集合J的所有可能子集的收集起來的集合

[desend_absU,index_desend] = sort(abs(U),'descend');%取絕對值後做降序排列
for ii = length(desend_absU):-1:1
    if desend_absU(ii)>1e-6%判斷productdes中非零值個數
        break;
    end
end
%Identify:Choose a set J of the K biggest coordinates
if ii>=Kin
    J = index_desend(1:Kin);%集合J (尚未正則化)
    Jval = desend_absU(1:Kin);%集合J對應的序列值，降冪排序後的值
    K = Kin;
else%or all of its nonzero coordinates,whichever is smaller
    J = index_desend(1:ii);%集合J(尚未正則化)
    Jval = desend_absU(1:ii);%集合J對應的序列值
    K = ii;
end
Emax = -1; %初始化能量值
for kk = 1 : K
    Et = Jval(kk)^2;
    mm = kk;
    Do = true;
    while Do
        mm = mm + 1;
        if mm > K
            mm = K;
            break;
        end
        Do = (Jval(kk) <= (2*Jval(mm)));
        if Do
            Et = Et + Jval(mm)^2;
        else
            mm = mm - 1;
        end
    end
    if Et > Emax
        Emax = Et;
        J0 = J(kk:mm);
    end
end

end
