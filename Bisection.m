function [rootBiSec] = Bisection(leftEnd, rightEnd, polyCof)
%% 函数解释：两分法求根
leftEndValue = PolynomialValue(polyCof, leftEnd);
rightEndValue = PolynomialValue(polyCof, rightEnd);

if abs(leftEndValue) < 10^-14
    rootBiSec = leftEnd;
    return;
end
if abs(rightEndValue) < 10^-14
    rootBiSec = rightEnd;
    return;
end
if leftEndValue*rightEndValue > 0
    display('Bisection failed, wrong input.');
    rootBiSec = 0;
    return;
end
mid = (leftEnd + rightEnd)/2;
midValue = PolynomialValue(polyCof, mid);

k = 0;

while abs(midValue) > 10^-12 && k < 500
    k = k + 1;
    if leftEndValue*midValue > 0
        leftEnd = mid;
    else
        rightEnd = mid;
    end
    
    leftEndValue = PolynomialValue(polyCof, leftEnd);
    rightEndValue = PolynomialValue(polyCof, rightEnd);
    mid = (leftEnd + rightEnd)/2;
    midValue = PolynomialValue(polyCof, mid);
end
rootBiSec = mid;