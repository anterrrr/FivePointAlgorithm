%% 本脚本的命名方式  变量：以小写字母开头，采用驼峰命名法  函数：以大写字母开头，通常为一个单词，当是多个单词是采用驼峰命名法  循环变量：以temp替代
clear;clc;close all;
%% 数据初始化
intrinsicParameters = load('./data/inParams.txt');camera1 = load('./data/codeCrd2D_Cam1.txt');camera2 = load('./data/codeCrd2D_Cam2.txt');
xp = intrinsicParameters(1);yp = intrinsicParameters(2);f = intrinsicParameters(3);
% 取前5对匹配点
cameraLeft = camera1(1:5,:);cameraRight = camera2(1:5,:); 
cameraLeftRest = camera1(6:end,:);cameraRightRest = camera2(6:end,:);

%% 畸变矫正
cameraLeft = cameraLeft-ones(length(cameraLeft),1)*[xp,yp]+Distortion(intrinsicParameters, cameraLeft);
cameraRight = cameraRight-ones(length(cameraRight),1)*[xp,yp]+Distortion(intrinsicParameters, cameraRight);
cameraLeftRest = cameraLeftRest-ones(length(cameraLeftRest),1)*[xp,yp]+Distortion(intrinsicParameters, cameraLeftRest);
cameraRightRest = cameraRightRest-ones(length(cameraRightRest),1)*[xp,yp]+Distortion(intrinsicParameters, cameraRightRest);

%% 归一化像面坐标
cameraLeft =  Normalization(intrinsicParameters,cameraLeft);
cameraRight =  Normalization(intrinsicParameters,cameraRight);
cameraLeftRest =  Normalization(intrinsicParameters,cameraLeftRest);
cameraRightRest =  Normalization(intrinsicParameters,cameraRightRest);

%% 求解基础解系
A = [];
[row,~] = size(cameraLeft);
for i = 1:row
    Cw2 = [cameraLeft(i,1)*cameraRight(i,1),cameraLeft(i,2)*cameraRight(i,1),cameraRight(i,1),cameraLeft(i,1)*cameraRight(i,2),cameraLeft(i,2)*cameraRight(i,2),cameraRight(i,2),cameraLeft(i,1),cameraLeft(i,2),1];
    A = [A;Cw2];
end

B = null(A);
Ew = B(:,1);
Ex = B(:,2);
Ey = B(:,3);
Ez = B(:,4);
Ew = reshape(Ew,[3,3]);
Ex = reshape(Ex,[3,3]);
Ey = reshape(Ey,[3,3]);
Ez = reshape(Ez,[3,3]);
Ew = Ew';
Ex = Ex';
Ey = Ey';
Ez = Ez';
syms w x y 
e = w*Ew+x*Ex+y*Ey+Ez;
e = vpa(e);

%% 求解系数矩阵cw
% 由det(E)=0得到第一个方程，进而得到1*10的系数矩阵
deltaE = det(e);
X = [y^3;x*y^2;x^2*y;x^3;y^2;x*y;x^2;y;x;1];
cw1 = [];
%将多项式分解成矩阵乘积的形式
for i=1:length(X)
    [r,q] = polynomialReduce(deltaE,X(i));
    deltaE = r;
    cw1 = [cw1,q];
end
%对多项式进行降幂排列
cw1 = Reorder(cw1, length(X)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %验证Cw1中的数据是否正确
% X = [1;x;y;x^2;x*y;y^2;x^3;x^2*y;x*y^2;y^3];
% sum1 = 0;
% for i=1:n
%     temp = Cw1(1,i)*X(i);
%     sum1 = sum1+temp;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 由2EE'E-trace(EE')E=0得到后九个方程--进而得到9*10的系数矩阵
D = 2*e*e.'*e-trace(e*e.')*e;
[row,col] = size(D);
cw2 = [];
for i=1:row
    for j=1:col
        temp = [];
        for k=1:length(X)
            [r,q] = polynomialReduce(D(i,j),X(k));
            D(i,j) = r;
            temp = [temp,q];
        end
        cw2 = [cw2;temp];
    end 
end


% % 对cw2中的每一行进行降幂排列
[row,col] = size(cw2);
cw2_ = [];
for i=1:row
    temp = cw2(i, :);
    temp = Reorder(temp, length(X));
    cw2_ = [cw2_; temp];
end
cw2 = cw2_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %验证cw2中的数据是否正确
% X = [1;x;y;x^2;x*y;y^2;x^3;x^2*y;x*y^2;y^3];
% sum2 = 0;
% for i=1:10
%     temp = cw2(8,i)*X(i);%要验证第几行则将相应的行数设置为几
%     sum2 = sum2+temp;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%将Cw1和Cw2合并
cw = [cw1;cw2];

%% 使用内置函数det计算系数矩阵cw的行列式
deltaCw = expand(det(cw));
disp('内置函数det计算|cw|结果为:');disp(deltaCw);

%% 使用Toeplitz矩阵计算系数矩阵cw的行列式
%将Cw中的每一个元素都用向量表示
V = cell(10);
X = [w^10;w^9;w^8;w^7;w^6;w^5;w^4;w^3;w^2;w;1];
for i=1:10
    for j=1:10
        temp = PolynomialVector(cw(i,j),X);
        V(i,j) = {temp};
    end
end

%计算每一个多项式的Toeplitz矩阵
T = cell(10);
for i=1:10
    for j=1:10
        temp = cell2mat(V(i,j));
        temp = ToeplitzMatrix(temp);
        T(i,j) = {temp};
    end
end
%使用Toeplitz矩阵运算求解Cw的行列式的表达式
n = 10;
G = cell(10);%G中存放用来描述多项式的矩阵
X = [w^10;w^9;w^8;w^7;w^6;w^5;w^4;w^3;w^2;w;1];
for k=1:n-1
    for i = k+1:n
        for j=k+1:n
            temp1=cell2mat(T(k,k));
            temp2=cell2mat(V(i,j));
            temp3=cell2mat(T(i,k));
            temp4=cell2mat(V(k,j));
            temp=(temp1*temp2-temp3*temp4);
            if k==1
                temp = temp(11:21,:);
            else      
                temp = inv(m.'*m)*m.'*temp;
            end
            temp5 = temp;
            temp6 = ToeplitzMatrix(temp5);
            V(i,j) = {temp5};
            T(i,j) = {temp6};
        end
    end
    m = cell2mat(T(k,k));
end
%计算得到Cw的行列式，行列式的向量表示即为V(10,10)
cwVector = cell2mat(V(10,10));
toeplitzCw = vpa(cell2mat(V(10,10)).'*X);
disp('Toeplitz矩阵计算|cw|结果为:');disp(toeplitzCw);

%% Sturm定理确定实根个数与单根区间
% 注：用Toeplitz矩阵运算得到系数矩阵cw的行列式 和 调用内置函数det所得结果几乎完全一致，多项式系数在小数点后面20位才有变化，
% 因此，笔者在本次算法复现中仍使用Toeplitz求解得到的行列式
%计算多项式|cw|的sturm序列
p0(w) = toeplitzCw;
p1(w) = diff(p0(w),w);
sturm = [];
X = [w^10;w^9;w^8;w^7;w^6;w^5;w^4;w^3;w^2;w;1];
count = 0;
i = 1;
while i==1
    [r,q] = polynomialReduce(p0,p1);
    r = -r;
    p0 = p1;
    p1 = r;
    if r==0
        break
    end
    temp = PolynomialVector(r,X);
    sturm = [sturm;temp.'];
    count = count+1;
end
p0(w) = toeplitzCw;
p1(w) = diff(p0(w),w);
sturm = [PolynomialVector(p0,X).';PolynomialVector(p1,X).';sturm];
%% testSturm+递归二分法求解多项式的根
a = -1000; b = 1000;
[rootIntervals] = RootIntervalsDetect(a, b, sturm);%单根区间
[numRoots] = NumRootsSturm(a, b, sturm);
realRoots = zeros(numRoots, 1);

for iRoots = 1 : numRoots
    [rootBiSec] = Bisection(rootIntervals(iRoots, 1), rootIntervals(iRoots, 2), toeplitzCw);
    realRoots(iRoots, 1) = rootBiSec;
end
disp('Sturm+递归二分法求解多项式实根为：');disp(realRoots);

%% roots函数求解多项式的根
%使用sturm序列确定多项式根的数目
X = [w^10;w^9;w^8;w^7;w^6;w^5;w^4;w^3;w^2;w;1];
sturmPolynomial(w) = vpa(sturm*X);
a = -1e+20;%分别使用a,b代表负无穷和正无穷
b = +1e+20;
a = sturmPolynomial(a);
b = sturmPolynomial(b);
ca = 0;
cb = 0;
[row,~] = size(sturm);

for i = 1:row-1
    if a(i)*a(i+1)<0
        ca = ca+1;
    end
    if b(i)*b(i+1)<0
        cb = cb+1;
    end
end
rootCounts1 = abs(ca-cb);

%matlab内置函数roots可以求解多项式的复根
complexRoots = roots(cell2mat(V(10,10)));
[row,~]=size(complexRoots);
for i=1:row
    complexRoots(i) = isreal(complexRoots(i));
end
rootCounts0 = sum(complexRoots);

disp(['内置函数roots解得多项式|Cw|实根个数为:',num2str(rootCounts0)])
disp(['sturm序列确定多项式|Cw|实根个数为:',num2str(rootCounts1)])
%在单根区间内使用二分法找到多项式|Cw|的数值解
disp('注：求解多项式实根时直接调用matlab内置函数roots，未采用文章中提到的二分法寻根')
disp('注：5点法真的难复现啊')

%matlab内置函数roots求解多项式实根
complexRoots = roots(cell2mat(V(10,10)));
realRoots = [];
[row,~]=size(complexRoots);
for i=1:row
    if double(isreal(complexRoots(i)))==1
        realRoots = [realRoots;complexRoots(i)];
    end
end
disp('内置函数roots求解多项式实根w为：');disp(realRoots);

%% 由w解算x,y
%将w带入cw中，可得到cw的数值形式
cwResult(w) = cw;
n = length(realRoots);
x = [];
y = [];
for i =1:n
    temp = cwResult(realRoots(i));
    temp_c = temp(:,2:10);%系数矩阵  
    temp_r = -temp(:,1);
    [U,S,V] = svd(temp_c);
    tempX = inv(V.')*pinv(S)*inv(U)*temp_r;
    x = [x;tempX(1)];
    x = double(x);
    y = [y;tempX(2)];
    y = double(y);
end
disp('解得x为：');disp(x);
disp('解得y为：');disp(y);

%% 计算得到所有E,并对E进行优化
e = [];
eAll = [];
errorAll = [];
[row,~] = size(cameraLeftRest);

for i=1:n
    tempE = realRoots(i)*Ew+x(i)*Ex+y(i)*Ey+Ez;
    %disp('E为：');disp(tempE);
    e = [e;tempE];
    [eOut,count] = OptimizationE(tempE,cameraLeftRest,cameraRightRest);
   error = [];
   
   %用优化后的E求解误差
   for j=1:length(cameraLeftRest)
       tempError = cameraRightRest(j,:)*eOut*cameraLeftRest(j,:).';
       error = [error;tempError];
   end
   
   eAll = [eAll;eOut];
   errorAll = [errorAll,error];
end

errorAll = abs(errorAll);
errorStd = std(errorAll);
errorMean = mean(errorAll);
errorMax = max(errorAll);



%由误差确定正确的本质矩阵E
index = [];
for i=1:n
    if errorMean(i)<1e-4
        index = [index,i];
    end
end

eTrue = [];
n = length(index);
for i=1:n
    temp = eAll(3*index(i)-2:3*index(i),:);
    eTrue = [eTrue;temp];
end



























