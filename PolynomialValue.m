function [y] = PolynomialValue(polyCof, v)
    syms w
    d(w) = polyCof;
    y = d(v);
end