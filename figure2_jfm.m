% Parameters
N = 100; % Number of time steps/planks
pi_m = 1; % Ratio of cylinder to liquid mass
pi_g = 1000; % Non-dimensional g - inversely proportional to nu^2
xi = 0; % Initial x-coordinate
yi = 0; % Initial y-coordinate
xxf = [2,5,10]; % Final x-coordinate
yyf = [-3,-1,-1]; % Final y-coordinate
bessel = csvread('bessel.csv'); % Roots of 1st order bessel function of first kind

% Initial guess
dti = 5; % Initial guess for time step
gi = zeros(1,N);  
for i = 1:N
	gi(i) = 0.5; % Initial guess for angles
end

for iloop = 1:3
    xf = xxf(iloop);
    yf = yyf(iloop);
    % Optimizer
    x0 = [gi,dti];
    A = []; 
    b = [];
    Aeq = [];
    beq = [];
    lb = -1*ones(1,N+1)*pi/2; % lower bound
    ub = ones(1,N+1)*pi/2; % upper bound
    lb(N+1) = 0;
    ub(N+1) = 1000000; % upper bound on time - very high value
    nonlcon = @(x)constraint(x,N,xi,yi,xf,yf,pi_m,pi_g,bessel); % non linear constraints
    options = optimoptions('fmincon','Display','iter','maxFunEvals',100000); 
    fun = @(x)N*x(N+1); % objective function
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    disp(newline);
    disp('Optimisation complete');

    % Post-processing
    omega = zeros(1,N+1); 
    % Calculating the angular velocity at grid points
    for i = 1:N
        nl1 = 0;
        for j = 1:i-1
            nl1 = nl1 + coeff(i,j,x(N+1),bessel)*(omega(j+1)-omega(j));           
        end
        omega(i+1) = (((-4*x(N+1)*nl1)+(1+pi_m)*pi_g*sin(x(i))*x(N+1))/(1+2*pi_m)) + omega(i);
    end
    disp(newline);
    disp('angular velocity calculated');
    
    xcoord = zeros(1,N+1); % x-coordinates of brachistrochrone curve
    ycoord = zeros(1,N+1); % y-coordinated of brachistrochrone curve
    for i = 1:N
        xcoord(i+1) = xcoord(i) + (omega(i)+omega(i+1))*cos(x(i))*0.5*x(N+1);
        ycoord(i+1) = ycoord(i) - (omega(i)+omega(i+1))*sin(x(i))*0.5*x(N+1);
    end    
    [xpoint,ypoint] = cycloid(xf,yf); % Brachistrochrone for a point mass - cycloid 
    plot(xcoord,ycoord,xpoint,ypoint,'--','Linewidth',1.5);
    ylabel({'$Y$'},'FontUnits','points','interpreter','latex','FontSize',22,'FontName','Times');
    xlabel('$X$','FontUnits','points','interpreter','latex','FontSize',22,'FontName','Times');
    pbaspect([1 1 1]);
    hold on;
end