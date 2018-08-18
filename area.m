function ar = area(xf,yf,xcoord,ycoord) 
    ar = 0;
    t = -2*pi+0.5;
    while 1
        a = yf*(t-sin(t))-xf*(1-cos(t));
        b = yf*(1-cos(t))-xf*(sin(t));
        t = t-(a/b);
        check = abs(a);
        if check < 0.00000001
            break
        end
    end 
    j = 1;
    c = yf/(1-cos(t));
    num_int = 10000;
    for i = 1:num_int
        t1 = (i-1)*t/num_int;
        t2 = i*t/num_int;
        x1 = c*(t1-sin(t1));
        x2 = c*(t2-sin(t2));
        y1 = c*(1-cos(t1));
        if (xcoord(j+1) < x1)
            j = j+1;
            disp(x1)
        end
        ytemp = ycoord(j)+ (ycoord(j+1)-ycoord(j))*(x1-xcoord(j))/(xcoord(j+1)-xcoord(j));
        ar = ar + abs(ytemp-y1)*(x2-x1);
    end
end