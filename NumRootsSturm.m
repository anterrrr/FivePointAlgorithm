function [numRooTs] =  NumRootsSturm(a,b,sturm)
    syms w
    X = [w^10;w^9;w^8;w^7;w^6;w^5;w^4;w^3;w^2;w;1];
    sturm_polynomial(w) = vpa(sturm*X);
    a = sturm_polynomial(a);
    b = sturm_polynomial(b);
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
   numRooTs = abs(ca-cb);
end