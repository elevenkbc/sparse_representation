% K mean example
% �s�y�@��2��observation�A���|�Ӹs�E�I�A�s�E�I�����߬�(1,2), (5,6), (8,8), (7, 1)
% observation ���ƥج� 4*N_group��
% solving  
% min_{D,X}|| Y - DX ||_F  subject to X_i = e_k  for some k.
% �`�N�@K-mean ���@�w�|���ġAK-mean�O�_���Ĩ��M���l��

close all
clear all
clc

%�سy�@��test data �䤤���I�b (1,2), (5,6), (8,8), (7,1).
Closure = [1,2; 5,6; 8,8; 7,1];
N_group = 100;
Y = zeros(2, 4*N_group);
%�N�Iø�s�b�ϧΤW

Color_map = [0,0.4,0.2;0,0,1;0.6,0,0.6;1,0,0.2];
for i = 1 : 4
    temp = [Closure(i,1)'*ones(1,N_group); Closure(i,2)'*ones(1,N_group)];
    Y(:,((i-1)*N_group + 1) : i*N_group) = temp + 0.6*randn(2,N_group);
end
sz = 5;%�C���I���b�|�Ѽ�
OBJ = scatter(Y(1, :), Y(2, :), sz, 'k', 'filled');%�N����Iø�s��ϧΤW



C0 = 10*rand(2, 4); %C0 ����ldictionary
%�s�@ø�s�ʵe������

hold on
K_center = scatter(C0(1,:), C0(2,:), 65, Color_map, '^');

pause(1)

for i = 1 : 20

    [C0, X] = one_step_K_mean(Y, C0);
    %���s�Ndata �̤����W��
    Color_j = [];
    for j = 1 : 4*N_group
        for k = 1 : 4
            if X(k,j) == 1
                Color_j = [Color_j; Color_map(k,:)];
            end
        end
    end
    set(OBJ, 'CData', Color_j);
    set(K_center, 'XData', C0(1,:), 'YData', C0(2,:));
    drawnow;

    pause(0.005)
    S = Y - C0*X;
    err = norm(S(:));
    xlabel(['iteration = ', num2str(i), '  ||Y-DX||_F^2 =' num2str(err)]);
    title('K-mean Animation');
    pause(0.1);
end
