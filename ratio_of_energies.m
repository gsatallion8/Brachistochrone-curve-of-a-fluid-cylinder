N = 200; % Number of time steps/planks
xi = 0; % Initial x-coordinate
yi = 0; % Initial y-coordinate
xf = 200; % Final x-coordinate
yf = -1; % Final y-coordinate
Nprec = 100;
bessel = csvread('bessel.csv'); % Roots of 1st order bessel function of first kind
lim = 2;
potential_energy = zeros(1,lim);
phi = zeros(1,lim);
fluid_translational = zeros(1,lim);
fluid_rotational = zeros(1,lim);
cylinder_kinetic = zeros(1,lim);
time = zeros(1,lim);
a = zeros(1,lim);
viscous_dissipation = zeros(1,lim);
xcoord = zeros(lim,N+1);
ycoord = zeros(lim,N+1);
omega = zeros(lim,N+1);
gamma = zeros(lim,N);
dti = 0.002747;
gi = atan(-yf/xf)*ones(1,N);
t = zeros(1,lim);
del = zeros(1,lim);
dim_time = zeros(1,lim);
for i2 = 1:lim 
    tic;
    xf = 120*i2 - 95;
    pi_g = 1000;
    pi_m = 0.01;
    % Optimizer
    x0 = [gi,dti];
    A = []; 
    b = [];
    Aeq = [];
    beq = [];
    lb = -1*ones(1,N+1)*(pi/2); % lower bound
    ub = ones(1,N+1)*(pi/2); % upper bound
    lb(N+1) = 0;
    ub(N+1) = 1000000; % upper bound on time - very high value
    nonlcon = @(x)constraint(x,N,xi,yi,xf,yf,pi_m,pi_g,bessel); % non linear constraints
    options = optimoptions('fmincon','FunctionTolerance',0.00000001,'OptimalityTolerance',0.00000001,'ConstraintTolerance',0.00000001,'maxFunEvals',1000000); 
    fun = @(x)N*x(N+1); % objective function
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    dti = x(N+1);
    for i = 1:N
        gi(i) = x(i); % Initial guess for next iteration
        gamma(i2,i) = x(i);
    end    
    % Calculating the angular velocity at grid points
    for i = 1:N
        nl1 = 0;
        for j = 1:i-1
            nl1 = nl1 + coeff(i,j,x(N+1),bessel)*(omega(i2,j+1)-omega(i2,j));           
        end
        omega(i2,i+1) = (((-4*x(N+1)*nl1)+(1+pi_m)*pi_g*sin(x(i))*x(N+1))/(1+2*pi_m)) + omega(i2,i);
    end    
    for i = 1:N
        xcoord(i2,i+1) = xcoord(i2,i) + (omega(i2,i)+omega(i2,i+1))*cos(x(i))*0.5*x(N+1);
        ycoord(i2,i+1) = ycoord(i2,i) - (omega(i2,i)+omega(i2,i+1))*sin(x(i))*0.5*x(N+1);
    end
    a(i2) = area(xf,yf,xcoord(i2,:),ycoord(i2,:));
    idiss = zeros(1,N+1);
    tanvel = zeros(N+1,Nprec+1); 
    disp(i2);
    for j1 = 1:N+1
        for i = 0:Nprec
            tanvel(j1,i+1) = velocity(j1,x(N+1),i/Nprec,bessel,omega(i2,:));
        end
    end 
    for i = 0:N
        for j1 = 1:Nprec-1
            idiss(i+1) = idiss(i+1) + j1*(tanvel(i+1,j1+2)-tanvel(i+1,j1+1)*(1+1/j1))^2;
        end  
        idiss(i+1) = 2*idiss(i+1);
    end
    potential_energy(i2) = (1+pi_m)*pi_g*(yi-yf);        
    for i = 1:N
        phi(i2) = phi(i2) + (idiss(i)+idiss(i+1))*0.5*x(N+1);
    end
    phi(i2) = phi(i2)/potential_energy(i2);
    viscous_dissipation(i2) = phi(i2)/pi_g;
    time(i2) = x(N+1);
    dim_time(i2) = time(i2)*sqrt(pi_g);
    fluid_translational(i2) = (0.5*omega(i2,N+1)*omega(i2,N+1))/potential_energy(i2);
    for i = 0:Nprec-1
        fluid_rotational(i2) = fluid_rotational(i2) + (i*tanvel(N+1,i+1)*tanvel(N+1,i+1))/(Nprec*Nprec);
    end
    fluid_rotational(i2) = fluid_rotational(i2)/potential_energy(i2);
    cylinder_kinetic(i2) = (pi_m*omega(i2,N+1)*omega(i2,N+1))/potential_energy(i2);  
    del(i2) = 1 - (fluid_translational(i2)+fluid_rotational(i2)+cylinder_kinetic(i2)+phi(i2));
    t(i2) = toc;
    disp('Time taken for the last iteration is');
    disp(t(i2));
end
csvwrite('potential_energy1_xfx_1.csv',potential_energy);
csvwrite('fluid_translational1_xfx_1.csv',fluid_translational);
csvwrite('fluid_rotational1_xfx_1.csv',fluid_rotational);  
csvwrite('cylinder_kinetic1_xfx_1.csv',cylinder_kinetic);
csvwrite('phi1_xfx_1.csv',phi);
csvwrite('time1_xfx_1.csv',time);
csvwrite('area1_xfx_1.csv',a);
csvwrite('xcoord1_xfx_1.csv',xcoord);
csvwrite('ycoord1_xfx_1.csv',ycoord);
csvwrite('gamma1_xfx_1.csv',gamma);
csvwrite('del1_xfx_1.csv',del);
csvwrite('viscous_dissipation1_xfx_1.csv',viscous_dissipation);
csvwrite('dim_time1_xfx_1.csv',dim_time);
















