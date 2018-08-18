% Parameters
N = 100; % Number of time steps/planks
xi = 0; % Initial x-coordinate
yi = 0; % Initial y-coordinate
xf = 10; % Final x-coordinate
yf = -1; % Final y-coordinate
bessel = csvread('bessel.csv'); % Roots of 1st order bessel function of first kind
dti = 5; % Initial guess for time step
gi = atan(-yf/xf)*ones(1,N); 
max = zeros(2,11);
for iloop = 1:11
    ipm = -1 + (iloop-1)*0.2;
    pi_m = power(10,ipm);
    inloop = 1;
    t1 = 0;
    t2 = 100000;
    t3 = 0;
    p1 = 0;
    p2 = 0;
    p3 = 0; 
    while 1
        ipg = 2.7 + (inloop-1)*0.1;
        pi_g = power(10,ipg);  
        % Optimizer
        y0 = [gi,dti];
        A = []; 
        b = [];
        Aeq = [];
        beq = [];
        lb = -1*ones(1,N+1)*pi/2; % lower bound
        ub = ones(1,N+1)*pi/2; % upper bound
        lb(N+1) = 0;
        ub(N+1) = 1000000; % upper bound on time - very high value
        nonlcon = @(y)constraint_particle(y,N,xi,yi,xf,yf,pi_g); % non linear constraints
        options = optimoptions('fmincon','Display','iter','maxFunEvals',100000); 
        fun = @(y)N*y(N+1); % objective function
        y = fmincon(fun,y0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        p1 = p2;
        p2 = p3;
        p3 = y(N+1);
        % Optimizer
        x0 = y;
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
        t1 = t2;
        t2 = t3;
        t3 = x(N+1);
        if (t1<t2) && (t3<t2)
            max(1,iloop) = power(10,(2.7+(inloop-2)*0.1));
            max(2,iloop) = t2/p2;
            break
        end
        inloop = inloop+1;
    end
end
dlmwrite('figure_3c.txt',max);