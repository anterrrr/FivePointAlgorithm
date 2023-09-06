function [y]=ToeplitzMatrix(x)
%将一个向量转换成toeplitz矩阵
n =length(x);
y = zeros(2*n-1,n);
for i=1:2*n-1
    for j=1:n
        for k =1:n
            if i-j==k-1
                y(i,j)=x(k);
            end
        end
    end
end

end