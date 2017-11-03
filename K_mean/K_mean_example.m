% K mean example
% �s�y�@��2��observation�A���T�ӥ|�Ӹs�E�I�A�s�E�I�����߬�(1,2), (5,6), (8,8), (7, 1)
% observation ���ƥج� 40��
% solving  
% min_{D,X}|| Y - DX ||_F  subject to X_i = e_k  for some k.
% �`�N�@K-mean ���@�w�|���ġAK-mean�O�_���Ĩ��M���l��

close all
clear all
clc
Closure = [1,2; 5,6; 8,8; 7,1];
Y = zeros(2, 40);
for i = 1 : 4
    for j = 1 : 10
        Y(:,(i-1)*10 + j) = Closure(i,:)' + 0.4*randn(2,1);
    end
end
plot(Y(1,:), Y(2,:), '.');
hold on
axis([0 ,13, 0, 12])


C0 = 8*rand(2, 4); %C0 ����ldictionary
Dic = plot(C0(1,:), C0(2,:), '*r'); %�s�@ø�s�ʵe������
legend('data', 'cluster centers');
for i = 1 : 20
    set(Dic, 'xdata', C0(1,:), 'ydata', C0(2,:));    
    drawnow;
    [C0, X] = K_mean_implement(Y, C0, 1);
    
    S = Y - C0*X;
    err = norm(S(:));
    xlabel(['iteration = ', num2str(i), '  ||Y-DX||_F^2 =' num2str(err)]);
    title('K-mean Animation');
    pause(0.1);
end
