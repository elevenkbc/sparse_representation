function [A] = ROMP(D, X, L)
% Regularized orthogonal matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% X ��l�T����@column�Ʀ����T���x�} X = [x_1, x_2, x_3, ..., x_p];
% D �T���r�� (DCT, harr wavelet ����)
% L a���D�s���Ӽ�
% output :
% a ��ܫY��

[n, p] = size(X);
k = size(D, 2);
A = zeros(k, p);

%���ND ���C�@��atom����
norm_D = zeros(size(D));
for i = 1 : k
    norm_D(:,i) = D(:,i)/norm(D(:,i));
end

for i = 1 : p
    %initailize
    R = X(:, i);
    atom_ind = [];%��l�U�ж��X���Ŷ��X
    for j = 1 : L
        g = abs(norm_D'*R); %�p�⤺�n�����
        if (sum(g~=0) <= L) %�p�Gg���D�s���p��L�A�h����������U��
            t1 = (1:size(g,1)).*(g~=0); 
            current_index = t1(t1~=0);
        else
            %��X���n����ȳ̤j��L�Ӥ��q�A�åB�NL�Ӥ��q�������U�Цs�� max_L_index
            [sorted_g , total_index]  =sort(g ,'descend');
            max_L_index = total_index(1:L);
            current_index = max_L_index;
        end
        
        R_current_index = Regular_index(current_index, g); %�N�U�Х��h�ơA�o�@�B�|�j�M�Ҧ���power set�A�GL���i�H�Ӥj
        
        if length([atom_ind, R_current_index'])>= 2*L %�p�G |atom_ind|>= 2L�A���X�j��
            break;
        else
            atom_ind = [atom_ind, R_current_index']; %��s�U�ж�
            D_atom_ind = norm_D(:,atom_ind); %��s�U�ж��ҹ������r��
            
            PinvD = pinv(D_atom_ind);
            temp_coe = PinvD*X(:,i); %�p��N X(:,j) ��v�� R(D_atom_ind) �Ŷ��ϱo|| X(:,i) - Phi*(temp_coe)|| �̤p���Y��
            P = D_atom_ind*PinvD; %P ����v�� R(D_atom_ind) �Ŷ��������v�x�}�A��v�x�}�������� P = D*inv(D'*D)*D'
            R = X(:,i) - P*X(:,i); %R ���ѤU���ݮt�V�q�AR�ݺ����P D_atom_ind �C��column�����۫���
        end
    end
    %�N�Y�ƶ�^ �Y�Ưx�}��
    A(atom_ind, i) = temp_coe;
end
end
function [J0]= Regular_index(J, U)
%J0 �����h�ƫ᪺J�y�Фl��
% ���h�ƪ��y�Фl�� J0 ���H�U�̨Τư��D����
% J0 = arg max { sum_{j in J0} |U_j|^2 }
% subject to  |u(i)| <= 2 |u(j)| for all i, j in J0,  where J0 in P(J)
% P(J) ���J��power set�A�����XJ���Ҧ��i��l���������_�Ӫ����X

%�䤤 J, U ���O��V�q
length_J = length(J);
B = ff2n(length_J); %B���C�@��row �� J���@�ؤl���i��
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