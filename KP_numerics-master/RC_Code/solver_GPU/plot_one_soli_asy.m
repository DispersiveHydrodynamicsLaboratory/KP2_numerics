

% 	sa = 1; q = 1;
%     a     = @(x,t) sa^2 .* ones(size(x));
%     ax    = @(x,t) zeros(size(x));
%     q     = @(x,t)  q   .* ones(size(x));
%     qx    = @(x,t) zeros(size(x));
%     uasoli = @(theta,x,t,a) a(x,t).*(sech(sqrt(a(x,t)/12).*theta)).^2;
%   dxuasoli = @(theta,dxtheta,x,t,a,ax) sech(sqrt(a(x,t)/12)*theta).^2.*...
%                (ax(x,t)+sqrt(a(x,t)/12).*tanh(sqrt(a(x,t)/12)*theta).*...
%                     dxtheta);
%     theta = @(x,y,t,a,q) (x + q(x,t).*y - (a(x,t)/3+q(x,t).^2) * t);
%   dxtheta = @(x,y,t,a,ax,q,qx) (-(x+(y-t).*q(x,t))*ax(x,t)+a(x,t).*...
%                            (-2+t*ax(x,t)+2*(t-y).*qx(x,t)));
%     uasy  = @(x,y,t)  uasoli(theta(x,y,t,a,q),x,t,a);
%   dxuasy  = @(x,y,t) dxuasoli(theta(x,y,t,a,q),...
%                             dxtheta(x,y,t,a,ax,q,qx),...
%                              x,t,a,ax);
%     u0    = @(x,y)    uasy(x,y,0);


	sa = 1; q = 1;
    a     = @(x,y,t) sa^2 .* ones(size(x));
    ax    = @(x,y,t) zeros(size(x));
    ay    = @(x,y,t) zeros(size(x));
    q     = @(x,y,t)  q   .* ones(size(x));
    qx    = @(x,y,t) zeros(size(x));
    qy    = @(x,y,t) zeros(size(x));
    uasoli = @(theta,x,y,t,a) a(x,y,t).*(sech(sqrt(a(x,y,t)/12).*theta)).^2;
  dxuasoli = @(theta,dxtheta,x,y,t,a,ax) sech(sqrt(a(x,y,t)/12).*theta).^2.*...
               (ax(x,y,t) + sqrt(a(x,y,t)/12).*tanh(sqrt(a(x,y,t)/12).*theta).*...
                    dxtheta);
  dyuasoli = @(theta,dytheta,x,y,t,a,q) sech(sqrt(a(x,y,t)/12).*theta).^2.*...
               (ay(x,y,t) + sqrt(a(x,y,t)/12).*tanh(sqrt(a(x,y,t)/12).*theta).*...
                    dytheta);

    theta = @(x,y,t,a,q) (x + q(x,y,t).*y - (a(x,y,t)/3+q(x,y,t).^2) * t);
  dxtheta = @(x,y,t,a,ax,q,qx) (-(x+(y-t.*q(x,y,t)).*q(x,y,t))*ax(x,y,t)+...
                    a(x,y,t).*(-2+t*ax(x,y,t)-2*(y-2*t*q(x,y,t)).*qx(x,y,t)));
  dytheta = @(x,y,t,a,ay,q,qy) (-(x+(y-t.*q(x,y,t)).*q(x,y,t)).*ay(x,y,t)+...
                    a(x,y,t).*(t.*ay(x,y,t)-2*y.*qy(x,y,t)+...
                    q(x,y,t)*(-2+4*t.*qy(x,y,t))));
  
                       
    uasy  = @(x,y,t)  uasoli(theta(x,y,t,a,q),x,y,t,a);
  dxuasy  = @(x,y,t) dxuasoli(theta(x,y,t,a,q),...
                            dxtheta(x,y,t,a,ax,q,qx),...
                             x,y,t,a,ax);
  dyuasy  = @(x,y,t) dyuasoli(theta(x,y,t,a,q),...
                            dxtheta(x,y,t,a,ax,q,qx),...
                              x,y,t,a,q);

    u0    = @(x,y)    uasy(x,y,0);
    
    

	  Lx = 50; Ly = Lx;
      Nx = 2^6; Ny = Nx;
      xplot  = (2*Lx/Nx)*[-Nx/2:Nx/2];
      yplot  = (2*Ly/Ny)*[-Ny/2:Ny/2];
      [XPLOT,YPLOT] = meshgrid(xplot,yplot);
      tplot = 0;
      u_init = u0(XPLOT,YPLOT);
      
        % Plot initial conditions and boundary conditions
        fontsize = 12;
        figure(1); clf;
        subplot(2,2,1)
            contourf(XPLOT,YPLOT,u_init,100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
            title('Initial Conditions');
        subplot(2,2,2)
            contourf(XPLOT,YPLOT,uasy(XPLOT,YPLOT,0),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
            title('Asymptotic u');
        subplot(2,2,3)
            contourf(XPLOT,YPLOT,dxuasy(XPLOT,YPLOT,0),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
            title('Asymptotic u, x-deriv');
        subplot(2,2,4)
            contourf(XPLOT,YPLOT,dyuasy(XPLOT,YPLOT,0),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
            title('Asymptotic u, y-deriv');
        set(gca,'fontsize',fontsize,'fontname','times');
        drawnow;
