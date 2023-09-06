function [EOut,n] = OptimizationE(E,cam1,cam2)
%% 函数解释对本质矩阵E进行局部优化
%E：本质矩阵E cam1：相机1归一化坐标 cam2:相机2归一化坐标 
%E_out:优化后的本质矩阵 n:迭代次数

%% 参数初始化
[row,~]=size(cam1);
epsilon0 = 1e-4;%设定收敛判断条件

%% 对E进行奇异值分解
[U,I,V]=svd(E);
scale = (I(1,1)+I(2,2))/2;
% scale = sqrt(scale);

%% 循环迭代
%注：5个角度初值为0，相当于乘
n=0;%迭代次数
i=1;
while i==1

    J = [];
    epsilon=[];
    %构造误差项和雅各比矩阵J
    for j=1:row
        x2 = cam2(j,:)*U;
        x1 = V.'*cam1(j,:).';
        
        temp_epsilon = scale*(x1(1)*x2(1)+x1(2)*x2(2));
        epsilon=[epsilon;temp_epsilon];
        %Right hand model 
        temp_j = [-x1(2)*x2(3), x1(1)*x2(3), x2(1)*x1(2)-x1(1)*x2(2), x2(1)*x1(3), -x2(2)*x1(3)];
  
        temp_j = scale*temp_j;
        J = [J;temp_j];
    end
    
    %计算得到角度修正量
    theta = -inv(J.'*J)*J.'*epsilon;
    u=theta(1);v=theta(2);w=theta(3);y=theta(4);x=theta(5);
    
    %Right hand model
    Ru = [1,0,0;0,cos(u),sin(u);0,-sin(u),cos(u)];
    Rv = [cos(v),0,-sin(v);0,1,0;sin(v),0,cos(v)];
    Rw = [cos(w),sin(w),0;-sin(w),cos(w),0;0,0,1];
    Ry = [cos(y),0,-sin(y);0,1,0;sin(y),0,cos(y)];
    Rx = [1,0,0;0,cos(x),sin(x);0,-sin(x),cos(x)];


    %更新U V，并使用新的U V计算误差
    U = U*Ru*Rv*Rw;
    V = V*Rx*Ry;
    epsilon_new = [];
    for k=1:row
        temp = cam2(k,:)*U*I*V.'*cam1(k,:).';
        epsilon_new = [epsilon_new;temp];
    end
    epsilon_new = abs(max(epsilon_new,[],'all'));

    
    %更新迭代次数
    n = n+1;
   
    %当新的误差满足条件时跳出while循环
    if epsilon_new<epsilon0
        break
    end
    
    if n>30
        break
    end
    
end

EOut = U*I*V.';
end
