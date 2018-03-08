%% Script for debugging asymptotic solution

%% Numerical Parameters
tmax   = 5;      % Solver will run from t=0 to t = tmax
numout = 10*tmax+1; % numout times will be saved (including ICs)
Lx     = 25;     % Solver will run on x \in [-Lx,Lx]
Ly     = 15;     % Solver will run on y \in [-Ly,Ly]
Nx     = 2^9;    % Number of Fourier modes in x-direction
Ny     = 2^8;    % Number of Fourier modes in y-direction

t      = linspace(0,tmax,numout);


    %% One-soliton
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
  dxtheta = @(x,y,t,a,ax,q,qx) (-(x+(y-t.*q(x,y,t)).*q(x,y,t)).*ax(x,y,t)+...
                    a(x,y,t).*(-2+t*ax(x,y,t)-2*(y-2*t*q(x,y,t)).*qx(x,y,t)));
  dytheta = @(x,y,t,a,ay,q,qy) (-(x+(y-t.*q(x,y,t)).*q(x,y,t)).*ay(x,y,t)+...
                    a(x,y,t).*(t.*ay(x,y,t)-2*y.*qy(x,y,t)+...
                    q(x,y,t).*(-2+4*t.*qy(x,y,t))));
  
                       
    uasy  = @(x,y,t)  uasoli(theta(x,y,t,a,q),x,y,t,a);
  dxuasy  = @(x,y,t) dxuasoli(theta(x,y,t,a,q),...
                            dxtheta(x,y,t,a,ax,q,qx),...
                             x,y,t,a,ax);
  dyuasy  = @(x,y,t) dyuasoli(theta(x,y,t,a,q),...
                            dxtheta(x,y,t,a,ax,q,qx),...
                              x,y,t,a,q);

    u0    = @(x,y)    uasy(x,y,0);

    
            % Load initial data
          xplot  = (2*Lx/Nx)*[-Nx/2:Nx/2];
          yplot  = (2*Ly/Ny)*[-Ny/2:Ny/2];
          [XPLOT,YPLOT] = meshgrid(xplot,yplot);
          tplot  = linspace(0,tmax,floor(tmax*10));
          u_init = u0(XPLOT,YPLOT);
        
        % Plot initial conditions and boundary conditions
        fontsize = 12;
        figure(1); clf;
        subplot(2,2,1)
            contourf(XPLOT,YPLOT,u_init,100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
            title('Initial Conditions');
        for ti = linspace(0,tmax,round(tmax))
            subplot(2,2,2)
                contourf(XPLOT,YPLOT,uasy(XPLOT,YPLOT,ti),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title(['Asymptotic u, t=',num2str(ti)]);
            subplot(2,2,3)
                contourf(XPLOT,YPLOT,dxuasy(XPLOT,YPLOT,ti),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Asymptotic u, x-deriv');
            subplot(2,2,4)
                contourf(XPLOT,YPLOT,dyuasy(XPLOT,YPLOT,ti),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Asymptotic u, y-deriv');
%             set(gca,'fontsize',fontsize,'fontname','times');
            drawnow;
            pause(0.25);
        end
        
%% Windowing function (from Kao 2010), rescaled, and derivs
n = 27; an = (1.111)^n*log(10);
W.o = exp( -an * abs(YPLOT/Ly).^n );
W.p = -an*n/Ly * W.o .* abs(YPLOT/Ly).^(n-1).*sign(YPLOT/Ly);
W.pp = an*n/Ly^2 * W.o .* abs(YPLOT/Ly).^(n-2) .* ...
           ( (-(n-1) + an*n*abs(YPLOT/Ly).^n).*sign(YPLOT/Ly).^2 );
figure(2); clf;
    subplot(3,1,1)
        plot(YPLOT(:,1),W.o(:,1))
    subplot(3,1,2)
        plot(YPLOT(:,1),W.p(:,1),...
             YPLOT(2:end,1),diff(W.o(:,1))./diff(YPLOT(:,1)));
    subplot(3,1,3)
        plot(YPLOT(:,1),W.pp(:,1),...
             YPLOT(2:end,1),diff(W.p(:,1))./diff(YPLOT(:,1)));

        