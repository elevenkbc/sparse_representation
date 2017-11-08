function [C, X] = one_step_K_mean(Y, C0)
%K mean algorithm
% min_{X,C} {|| Y - CX||_F^2} subject to for all i x_i = e_k for some k
% iterN �|�N����
% ob_dim ���[��Ȫ����סA�]�N�OY��Row���ƥءA�`�NY�����C�@��column���@���[���
% N ���[��Ȫ��Ӽ�
% K ���r�� atom ���Ӽ�
% C0 ���|�N����l�r��
% X �䤤�C��cloumn �� e_k for some k

N = size(Y, 2);
K = size(C0, 2);
X = zeros(K, N); %��l��codebook�@�Y�Ưx�}
C = C0;

% sparse coding
% �N observation �̷� codebook �������A�`�@����K��
% �Y observation  || Y_j - C_k || < || Y_j - C_l ||  which l ~= k
% �h�N observation Y_j ���b C_k ���O��
R = cell(1, K); %��cell �ӼаO���O   R{1} ���@��column �V�q�A�V�q���C�@�Ӥ����� Y��column ���аO
%�|�Ҩӻ��A R{1} = [1, 5, 6]  �Х� {Y_1, Y_5, Y_6} �ݩ� ����a��  C_1 �������C

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
%��������

% code book update
%�β�k�դ��� R_k �� Y �������ƨӧ�s�r�� �x�}C
X = zeros(K, N); %��l��codebook�@�Y�Ưx�}
for k = 1 : K
    R_k = R{k};
    if isempty(R_k)
        continue; %�Y�Ӥ����O�Ū��A�h���L
    end
    temp_sum = 0;
    for j = 1 : length(R_k)
        
        X(k, R_k(j)) = 1; %�p��X
        temp_sum = temp_sum + Y(:, R_k(j));
        
    end
    
    C(:,k) = temp_sum/length(R_k);
    
end



end