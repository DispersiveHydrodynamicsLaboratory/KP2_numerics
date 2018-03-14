t = 75; ti = 1;    
figure(3); clf;
    a = @(x,y,t) 9/64*(2-(y)/t).^2 .*(-2/3*t<(y)).*((y)< 2*t) + ...9/16*(1+(y+50)/(2*t)).^2 .*( 2/3*t>y+50).*(y+50>-2*t) + ...
                 1                .*(-2/3*t>=y);
             
    q = @(x,y,t) +1                 .*( 2*t<=y)+ ...
                 (1-3/8*(2-(y)/t)).*(-2/3*t<y).*(y< 2*t);
%     ti=75
    U = a(X,Y,t(ti)).*sech( sqrt(a(X,Y,t(ti))/12).*...
         (X + q(X,Y,t(ti)).*Y - (a(X,Y,t(ti))/3 + q(X,Y,t(ti)).^2)*t(ti))).^2;
    contourf(X,Y,U,50,'edgecolor','none'); 
    colormap(load('CoolWarmFloat257.csv')); colorbar;
