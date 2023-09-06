function [y] = PolynomialVector(x,X)
n = length(X);
temp = [];
for i=1:n
    [r,q] = polynomialReduce(x,X(i));
    x = r;
    q = double(q);
    temp = [temp;q];
end
y = temp;
end