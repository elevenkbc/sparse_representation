function [A] = ROMP(D, X, L)
% Regularized orthogonal matching pursuit  (greed algorithm)   
%  
%   min_{a}   || x - D a ||_2         subject to  ||a||_0  <= L
% 
%
% input :
% X ��l�T����@column�Ʀ����T���x�} X = [x_1, x_2, x_3, ..., x_p];
% D �T���r�� (DCT, harr wavelet ����) �C��column�������׬�1
% L a���D�s���Ӽ�
% output :
% a ��ܫY��

[n, p] = size(X);
k = size(D, 2);
A = zeros(k, p);


for i = 1 : p
    %initailize
    R = X(:, i);
    atom_ind = [];%��l�U�ж��X���Ŷ��X
    for j = 1 : L    
        R_current_index = Regular_index(D'*R, L); %�p�⥿�h�ƪ��y��
        
        if length([atom_ind, R_current_index'])>= 2*L %�p�G |atom_ind|>= 2L�A���X�j��
            break;
        else
            atom_ind = [atom_ind, R_current_index']; %��s�U�ж�
            D_atom_ind = D(:,atom_ind); %��s�U�ж��ҹ������r��
            
            PinvD = pinv(D_atom_ind);
            temp_coe = PinvD*X(:,i); %�p��N X(:,j) ��v�� R(D_atom_ind) �Ŷ��ϱo|| X(:,i) - Phi*(temp_coe)|| �̤p���Y��
            P = D_atom_ind*PinvD; %P ����v�� R(D_atom_ind) �Ŷ��������v�x�}�A��v�x�}�������� P = D*inv(D'*D)*D'
            R = X(:,i) - P*X(:,i); %R ���ѤU���ݮt�V�q�AR�ݺ����P D_atom_ind �C��column�����۫���
        end
        if norm(R) < 1.0e-6 %�p�G�ѤU���ݮt�q�Ӥp�A�N���X�j��
            break;
        end
    end
    %�N�Y�ƶ�^ �Y�Ưx�}��
    A(atom_ind, i) = temp_coe;
end
end
function [J0]= Regular_index(U, Kin)
%J0 �����h�ƫ᪺J�y�Фl��
% ���h�ƪ��y�Фl�� J0 ���H�U�̨Τư��D����
% J0 = arg max { sum_{j in J0} |U_j|^2 }
% subject to  |u(i)| <= 2 |u(j)| for all i, j in J0,  where J0 in P(J)
% P(J) ���J��power set�A�����XJ���Ҧ��i��l���������_�Ӫ����X

[desend_absU,index_desend] = sort(abs(U),'descend');%������ȫᰵ���ǱƦC
for ii = length(desend_absU):-1:1
    if desend_absU(ii)>1e-6%�P�_productdes���D�s�ȭӼ�
        break;
    end
end
%Identify:Choose a set J of the K biggest coordinates
if ii>=Kin
    J = index_desend(1:Kin);%���XJ (�|�����h��)
    Jval = desend_absU(1:Kin);%���XJ�������ǦC�ȡA�����Ƨǫ᪺��
    K = Kin;
else%or all of its nonzero coordinates,whichever is smaller
    J = index_desend(1:ii);%���XJ(�|�����h��)
    Jval = desend_absU(1:ii);%���XJ�������ǦC��
    K = ii;
end
Emax = -1; %��l�Ư�q��
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
