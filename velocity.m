function y = velocity(i,k,l,bessel,a)
y = 0;
for j = 1:i-1
    y1 = 1;
    y2 = 0;
    a1 = 1;
    while abs(y1) > power(10,-1)
        y1 = 2*exp(-1*bessel(a1)*bessel(a1)*(i-j)*k)*besselj(1,l*bessel(a1))/(bessel(a1)*besselj(0,bessel(a1)));
        y2 = y2 + y1;
        a1 = a1+1;
    end
    y = y + y2*(a(j+1)-a(j));
end
y = y +l*a(i);
end

