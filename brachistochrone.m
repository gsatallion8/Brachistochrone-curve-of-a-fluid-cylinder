pi_g=1000;
pi_m=0.01;
xi = 0; % Initial x-coordinate
yi = 0; % Initial y-coordinate
bessel = csvread('bessel.csv'); % Roots of 1st order bessel function of first kind
N=400;
xf = 100;
yf = -1;
% Initial guess
dti = 0.0002747; % Initial guess for time step
gi = atan(-yf/xf)*ones(1,N);
% Optimizer
x0 = [gi,dti];
% x0 = x;
A = []; 
b = [];
Aeq = [];
beq = [];
lb = -1*ones(1,N+1)*pi/2; % lower bound
ub = ones(1,N+1)*pi/2; % upper bound
lb(N+1) = 0;
ub(N+1) = 1000000000; % upper bound on time - very high value
nonlcon = @(x)constraint(x,N,xi,yi,xf,yf,pi_m,pi_g,bessel); % non linear constraints
options = optimoptions('fmincon','Display','iter','maxIter',100000,'maxFunEvals',10000000,'FunctionTolerance',0.0001,'OptimalityTolerance',0.00001,'ConstraintTolerance',0.00000001); 
fun = @(x)N*x(N+1); % objective function
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% Post-processing
omega = zeros(1,N+1); 
x1=zeros(1,N+1);
y1=zeros(1,N+1);
% Calculating the angular velocity at grid points
str=sprintf('Omega_%f_%f_%f_%f.txt',pi_m,pi_g,xf,yf);
fileID=fopen(str,'w');
fprintf(fileID,'%d %d\n',N,x(N+1)*N);
 for i = 1:N
    fprintf(fileID,'%f\n',omega(i));
   nl1 = 0;
    for j = 1:i-1
     nl1 = nl1 + coeff(i,j,x(N+1),bessel)*(omega(j+1)-omega(j));           
    end
    omega(i+1) = (((-4*x(N+1)*nl1)+(1+pi_m)*pi_g*sin(x(i))*x(N+1))/(1+2*pi_m)) + omega(i);
    x1(i+1)=x1(i)+(omega(i)+omega(i+1))*cos(x(i))*x(N+1)/2;
    y1(i+1)=y1(i)+(omega(i+1)+omega(i))*sin(x(i))*x(N+1)/2;
 end
fclose(fileID);

str=sprintf('Coordinates_%f_%f_%f_%f.txt',pi_m,pi_g,xf,yf);
fileID=fopen(str,'w');
for i = 1:N
    fprintf(fileID,'%f\t%f\n',x1(i), -y1(i));
end
fprintf(fileID,'%f\t%f\n',xf, yf);
fclose(fileID);

figure1=figure;
options = optimoptions('fsolve','Display','iter','maxIter',100000,'maxFunEvals',10000000,'FunctionTolerance',0.0001,'OptimalityTolerance',0.00001);
tf = fsolve(@(tf)cycloid_terminal(tf,xf,yf), 2*pi-0.02*pi, options);
r = xf/(tf - sin(tf));

t = linspace(0,tf,500);
xc = r*(t-sin(t));
yc = -r*(1-cos(t));

str=sprintf('Cycloid_%f_%f.txt',xf,yf);
fileID=fopen(str,'w');

plot(x1,-y1)
hold on
plot(xc,yc)

for i = 1:500
    if yc(i)<=yf && xc(i)<=xf
        fprintf(fileID,'%f\t%f\n',xc(i), yc(i));
    end
end

fprintf(fileID,'%f\t%f\n',xf, yf);
fclose(fileID);

string=sprintf('Brach_%d_%d.png',xf,yf);
saveas(figure1,string)