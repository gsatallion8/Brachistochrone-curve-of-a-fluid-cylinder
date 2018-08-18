FileName = 'Brachistochrone.gif';

xf = 100;
yf = -1;

for k = 1:2:401
    
    str = sprintf('Image_%f_%f_%d.jpg',xf,yf,k);
    img = imread(str);
    [imind,cm] = rgb2ind(img,256);
    if k ==1
        imwrite(imind,cm,FileName,'gif','LoopCount',Inf,'DelayTime',0.05);
    else
        imwrite(imind,cm,FileName,'gif','WriteMode','append','DelayTime',0.05);
    end
end