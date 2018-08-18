pi_g=1000;
pi_m=0.01;

xf = 100;
yf = -1;

options = optimset('Display','iter', 'TolX', 0.001);
fun = @(tf)cycloid_terminal(tf,xf,yf);
tf = fsolve(fun, 2*pi-0.02*pi, options);
r = xf/(tf - sin(tf));

str=sprintf('Coordinates_%f_%f_%f_%f.txt',pi_m,pi_g,xf,yf);

fileID=fopen(str,'r');
A = fscanf(fileID,'%f\t%f');
x = A(1:2:end);
y = A(2:2:end);

fclose(fileID);

yc = zeros(1,length(x));
t = pi;
for i = 1:length(x)
    fun = @(t)cycloid_interim(t,x(i),r);
    t = fsolve(fun, t, options);
    yc(i) = -r*(1-cos(t));
end

str=sprintf('Cycloid_%f_%f.txt',xf,yf);
fileID=fopen(str,'w');

for i = 1:length(x)
    fprintf(fileID,'%f\t%f\n',x(i), yc(i));
end

fclose(fileID);
    