function y = coeff(i,j,k,bessel) %sum of the infinte series 
    a = 1;
    y1 = 1;
    y = 0;
    while abs(y1) > power(10,-10)
        y1 = exp(-1*bessel(a)*bessel(a)*(i-j)*k);
        y = y + y1;
        a = a+1;
    end
end