pi_g=1000;
pi_m=0.01;
xf = 100;
yf = -1;

str=sprintf('Coordinates_%f_%f_%f_%f.txt',pi_m,pi_g,xf,yf);

fileID=fopen(str,'r');
A = fscanf(fileID,'%f\t%f');
x = A(1:2:end);
y = A(2:2:end);

fclose(fileID);

str=sprintf('Cycloid_%f_%f.txt',xf,yf);

fileID=fopen(str,'r');
A = fscanf(fileID,'%f\t%f');
xc = A(1:2:end);
yc = A(2:2:end);

fclose(fileID);

t = linspace(0,2*pi,400);

for i = 351:2:401
    
    disp(i);
    
    string = sprintf('Image_%f_%f_%d.jpg',xf,yf,i);
    figure1=figure('visible','off');
    plot(x, y, 'r');
    hold on
    plot(xc, yc, 'b');
    hold on
    
    X = zeros(2*i);
    Y = zeros(2*i);
    
    k=1;
    for j = 1:i
        X(k) = x(j);
        Y(k) = y(j);
        k=k+1;
    end
    for j = 1:i
        X(k) = xc(i+1-j);
        Y(k) = yc(i+1-j);
        k=k+1;
    end
    
    fill(X,Y,'g');
    hold on
    
    x_cir = x(i) + cos(t);
    y_cir = y(i) + sin(t);
    
    plot(x_cir, y_cir, 'k');
    
    axis([-10 110 -33.34 5])
    pbaspect([3 1 1])
    
    saveas(figure1,string);
end