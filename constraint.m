function [c,ceq] = constraint(x,N,xi,yi,xf,yf,pi_m,pi_g,bessel)
    o = zeros(1,N+1);
    for i = 1:N
       nl1 = 0;
       for j = 1:i-1
           nl1 = nl1 + coeff(i,j,x(N+1),bessel)*(o(j+1)-o(j));           
       end
       o(i+1) = (((-4*x(N+1)*nl1)+(1+pi_m)*pi_g*sin(x(i))*x(N+1))/(1+2*pi_m)) + o(i);
    end    
    c = [];
    ceq = zeros(1,2);
    nl2 = 0;
    nl3 = 0;
    for i = 1:N
       nl2 = nl2 + (o(i)+o(i+1))*0.5*cos(x(i));
       nl3 = nl3 + (o(i)+o(i+1))*0.5*sin(x(i));
    end   
    ceq(1) = xf-xi-nl2*x(N+1);
    ceq(2) = yf-yi+nl3*x(N+1); 
end