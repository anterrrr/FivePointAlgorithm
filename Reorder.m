function [reorder]=Reorder(cw, n)
%% 函数解释：对多项式进行降幂排列
for i=1:n
    temp = cw(1,i);
    cw(1,i) = cw(1,n+1-i);
    cw(1,n+1-i) = temp;
    if i==n/2
        break
    end
end
reorder = cw;
end