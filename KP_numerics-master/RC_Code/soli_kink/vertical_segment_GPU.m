function [ soli ] = vertical_segment_GPU( sau, sad, qau, qad, x0, y0, Lx, w )
%Creates a structure that's a soliton IC
% Inputs
% sa: square root of the amplitude
% qa: slope of the angle of propagation

verbose = 0; % nonzero means plotting, text output, etc
ad = sad^2;
au = sau^2;
%% Define Soliton parameters and derivatives in x and y
    soli.tw    = 1/10; %tanh width factor
    soli.w     = w;
    soli.a     = @(x,y,t) gpuArray((ad-au)/2*tanh(soli.tw*(y-soli.w/2)) - (ad-au)/2*tanh(soli.tw*(y+soli.w/2)));
    soli.ax    = @(x,y,t) gpuArray(zeros(size(x)));
    soli.ay    = @(x,y,t) gpuArray((ad-au)/2*soli.tw*sech(soli.tw*(y-soli.w/2)).^2 - (ad-au)/2*soli.tw*sech(soli.tw*(y+soli.w/2)).^2);
    soli.q     = @(x,y,t) gpuArray((qad-qau)/2*tanh(soli.tw*(y-soli.w/2)) - (qad-qau)/2*tanh(soli.tw*(y+soli.w/2)));
    soli.qx    = @(x,y,t) gpuArray(zeros(size(x)));
    soli.qy    = @(x,y,t) gpuArray((qad-qau)/2*soli.tw*sech(soli.tw*(y-soli.w/2)).^2 - (qad-qau)/2*soli.tw*sech(soli.tw*(y+soli.w/2)).^2);
    soli.x0    = gpuArray(x0);
    soli.y0    = gpuArray(y0);
    soli.G = gpuArray(0); % background
    
%% Define Soliton function, derivatives
    soli.u  = @(theta,x,y,t,a,g) g + a(x,y,t).*(sech(sqrt(a(x,y,t)/12).*theta)).^2; % CORRECT
    soli.dx = @(theta,dxtheta,x,y,t,a,ax) sech(sqrt(a(x,y,t)/12).*theta).^2.*...
               (ax(x,y,t) + sqrt(a(x,y,t)/12).*tanh(sqrt(a(x,y,t)/12).*theta).*...
                    dxtheta);
    soli.dy = @(theta,dytheta,x,y,t,a,q,ay) sech(sqrt(a(x,y,t)/12).*theta).^2.*...
               (ay(x,y,t) + sqrt(a(x,y,t)/12).*tanh(sqrt(a(x,y,t)/12).*theta).*...
                    dytheta);
% Theta and derivatives
    soli.th = @(x,y,t,a,q,g) (x + q(x,y,t).*y - (a(x,y,t)/3+q(x,y,t).^2+g) * t); % FIXED
    soli.thx = @(x,y,t,a,ax,q,qx,g) ((g*t-x-y.*q(x,y,t)+t*q(x,y,t).^2).*ax(x,y,t)+...
                        a(x,y,t).*(-2+t*ax(x,y,t)-2*(y-2*t*q(x,y,t)).*qx(x,y,t)));
    soli.thy = @(x,y,t,a,ay,q,qy,g) ((g*t-x-y.*q(x,y,t)+t*q(x,y,t).^2).*ay(x,y,t)+...
                        a(x,y,t).*(t.*ay(x,y,t)-2*y.*qy(x,y,t)+...
                        q(x,y,t).*(-2+ 4*t.*qy(x,y,t))));
  
% Asymptotic soliton approximation
    soli.ua  = @(x,y,t)  soli.u(soli.th(x-soli.x0,y-soli.y0,t,soli.a,soli.q,soli.G),...
                                x-soli.x0,y-soli.y0,t,soli.a,soli.G);
    soli.uax  = @(x,y,t) soli.dx(soli.th(x-soli.x0,y-soli.y0,t,soli.a,soli.q,soli.G),...
                            soli.thx(x-soli.x0,y-soli.y0,t,soli.a,soli.ax,soli.q,soli.qx,soli.G),...
                             x-soli.x0,y-soli.y0,t,soli.a,soli.ax);
    soli.uay  = @(x,y,t) soli.dy(soli.th(x-soli.x0,y-soli.y0,t,soli.a,soli.q,soli.G),...
                                 soli.thy(x-soli.x0,y-soli.y0,t,soli.a,soli.ay,soli.q,soli.qy,soli.G),...
                                  x-soli.x0,y-soli.y0,t,soli.a,soli.q,soli.ay);
% Initial condition
    soli.u0    = @(x,y)    soli.ua(x,y,0);

        % Figure for debugging
        if verbose
            xplot = -Lx:0.5:Lx;
            yplot = -Lx/2:0.5:Lx/2;
            [XPLOT,YPLOT] = meshgrid(xplot,yplot);
            figure(1); clf;
            subplot(2,2,1)
                contourf(XPLOT,YPLOT,soli.u0(XPLOT,YPLOT),100,'edgecolor','none'); xlabel('x'); ylabel('y');
                title('Initial Conditions');
            for ti = linspace(0,5,2)
                subplot(2,2,2)
                    contourf(XPLOT,YPLOT,soli.ua(XPLOT,YPLOT,ti),100,'edgecolor','none'); xlabel('x'); ylabel('y');
                    title(['Asymptotic u, t=',num2str(ti)]);
                subplot(2,2,3)
                    contourf(XPLOT,YPLOT,soli.uax(XPLOT,YPLOT,ti),100,'edgecolor','none'); xlabel('x'); ylabel('y');
                    title('Asymptotic u, x-deriv');
                subplot(2,2,4)
                    contourf(XPLOT,YPLOT,soli.uay(XPLOT,YPLOT,ti),100,'edgecolor','none'); xlabel('x'); ylabel('y');
                    title('Asymptotic u, y-deriv');
    %             set(gca,'fontsize',fontsize,'fontname','times');
                drawnow;
                pause(0.25);
            end
        end

end

