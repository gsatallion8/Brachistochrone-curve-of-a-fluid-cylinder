function z = cycloid_terminal( tf ,xf, yf )
z = (yf*(tf - sin(tf)) - xf*(1-cos(tf)))^2;
end