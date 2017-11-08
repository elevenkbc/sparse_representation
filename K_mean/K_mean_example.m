% K mean example
% 製造一組2維observation，有四個群聚點，群聚點的中心為(1,2), (5,6), (8,8), (7, 1)
% observation 的數目為 4*N_group個
% solving  
% min_{D,X}|| Y - DX ||_F  subject to X_i = e_k  for some k.
% 注意　K-mean 不一定會收斂，K-mean是否收斂取決於初始值

close all
clear all
clc

%建造一個test data 其中心點在 (1,2), (5,6), (8,8), (7,1).
Closure = [1,2; 5,6; 8,8; 7,1];
N_group = 100;
Y = zeros(2, 4*N_group);
%將點繪製在圖形上

Color_map = [0,0.4,0.2;0,0,1;0.6,0,0.6;1,0,0.2];
for i = 1 : 4
    temp = [Closure(i,1)'*ones(1,N_group); Closure(i,2)'*ones(1,N_group)];
    Y(:,((i-1)*N_group + 1) : i*N_group) = temp + 0.6*randn(2,N_group);
end
sz = 5;%每個點的半徑參數
OBJ = scatter(Y(1, :), Y(2, :), sz, 'k', 'filled');%將資料點繪製到圖形上



C0 = 10*rand(2, 4); %C0 為初始dictionary
%製作繪製動畫的物件

hold on
K_center = scatter(C0(1,:), C0(2,:), 65, Color_map, '^');

pause(1)

for i = 1 : 20

    [C0, X] = one_step_K_mean(Y, C0);
    %重新將data 依分類上色
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
