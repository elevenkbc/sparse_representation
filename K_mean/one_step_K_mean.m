function [C, X] = one_step_K_mean(Y, C0)
%K mean algorithm
% min_{X,C} {|| Y - CX||_F^2} subject to for all i x_i = e_k for some k
% iterN 疊代次數
% ob_dim 為觀察值的維度，也就是Y的Row的數目，注意Y中的每一個column為一個觀察值
% N 為觀察值的個數
% K 為字典 atom 的個數
% C0 為疊代的初始字典
% X 其中每個cloumn 為 e_k for some k

N = size(Y, 2);
K = size(C0, 2);
X = zeros(K, N); %初始化codebook　係數矩陣
C = C0;

% sparse coding
% 將 observation 依照 codebook 做分類，總共分成K類
% 若 observation  || Y_j - C_k || < || Y_j - C_l ||  which l ~= k
% 則將 observation Y_j 分在 C_k 類別中
R = cell(1, K); %用cell 來標記類別   R{1} 為一個column 向量，向量中每一個元素為 Y的column 的標記
%舉例來說， R{1} = [1, 5, 6]  標示 {Y_1, Y_5, Y_6} 屬於 比較靠近  C_1 的分類。

for i = 1 : N
    dis_2norm =  1000;
    temp = 0;
    for k = 1 : K
        if norm(Y(:,i) - C(:,k)) < dis_2norm
            dis_2norm = norm(Y(:,i) - C(:,k));
            temp = k;
        end
    end
    R{temp} = [R{temp}; i];
end
%分類結束

% code book update
%用第k組分類 R_k 中 Y 的平均數來更新字典 矩陣C
X = zeros(K, N); %初始化codebook　係數矩陣
for k = 1 : K
    R_k = R{k};
    if isempty(R_k)
        continue; %若該分類是空的，則跳過
    end
    temp_sum = 0;
    for j = 1 : length(R_k)
        
        X(k, R_k(j)) = 1; %計算X
        temp_sum = temp_sum + Y(:, R_k(j));
        
    end
    
    C(:,k) = temp_sum/length(R_k);
    
end



end