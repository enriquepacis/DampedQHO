function hc=DrawCircle(haxis,x,y,R,LineColor,LineWidth,FillColor)

N=40;
theta=linspace(0,2*pi,N);
X=x+R*cos(theta);
Y=y+R*sin(theta);

hc=patch(X,Y,FillColor);


set(hc, 'Parent', haxis);
set(hc, 'EdgeColor', LineColor);
set(hc, 'LineWidth', LineWidth);
