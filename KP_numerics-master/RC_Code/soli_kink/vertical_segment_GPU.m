function [ soli ] = vertical_segment_GPU(am,ad,au,...
                                     qm,qu,qd,...
                                     x0,y0,w,Lx,Ly,Nx,Ny,tstart)
	soli.x  = linspace(-Lx,Lx-(2*Lx/Nx),Nx);
	soli.y  = linspace(-Ly,Ly-(2*Ly/Ny),Ny);
    [X,Y] = meshgrid(soli.x,soli.y);
    soli.Y = Y;
    soli.X = X;
	%% CURRENTLY ASSUMES am=1
    if am~=1
        disp('Setting middle state to 1...');
        am = 1;
    end
          %% Define Soliton parameters and derivatives in x and y
            soli.w     = (w);
            soli.a     = @(x,y,t) ones(size(x)).*(...
                                  au                         .*( (y-soli.w/2)./t >= 2-4/3*sqrt(au) ) + ...
                                  9/64*(2-(y-soli.w/2)./t).^2 .*(-2/3<(y-soli.w/2)./t).*((y-soli.w/2)./t<2-4/3*sqrt(au)) +...
                                  1                          .*((y-soli.w/2)./t<=-2/3).*((y+soli.w/2)./t>=+2/3) +...
                                  9/64*(2+(y+soli.w/2)./t).^2 .*(+2/3>(y+soli.w/2)./t).*((y+soli.w/2)./t>4/3*sqrt(au)-2) +...
                                  ad                         .*( (y+soli.w/2)./t <= 4/3*sqrt(ad)-2)...
                                  );
            soli.ax    = @(x,y,t) (zeros(size(x)));
            soli.ay    = @(x,y,t) zeros(size(x));
%                                   9/32*(2-(y-soli.w/2)./t)/t .*(-2/3<(y-soli.w/2)./t).*((y-soli.w/2)./t<2-4/3*sqrt(au)) +...
%                                   9/32*(2+(y+soli.w/2)./t)/t .*(+2/3>(y+soli.w/2)./t).*((y+soli.w/2)./t>4/3*sqrt(au)-2);
            soli.q     = @(x,y,t) ((1-sqrt(soli.a(x,y,t)))     .*(-2/3<(y-soli.w/2)./t) + ...
                                  (sqrt(soli.a(x,y,t))-1)    .*(+2/3>(y+soli.w/2)./t));
            soli.qx    = @(x,y,t) (zeros(size(x)));
            soli.qy    = @(x,y,t) zeros(size(x));
            soli.x0    = (x0);
            soli.y0    = (y0);
            soli.G = 0;
      %% Define Soliton function, derivatives
          soli.u  = @(theta,x,y,t,a,g) g + a(x,y,t).*(sech(sqrt(a(x,y,t)/12).*theta)).^2; % CORRECT
          soli.dx = @(x,y,t) (zeros(size(x)));
          soli.dy = @(x,y,t) (zeros(size(x)));

          soli.th.intx = soli.X;
            
          soli.th.inty = zeros(Ny,Nx);
            for Nyi = 1:length(soli.y)
                soli.th.inty(Nyi,:) = integral(@(y)soli.q(0,y,tstart),0,soli.y(Nyi));
            end
            
          soli.t = 0.01:0.01:tstart;
          soli.th.intt = zeros(1,length(soli.t));
            for Nti = 1:length(soli.t)
                soli.th.intt(Nti) = -integral(@(t)(soli.a(0,0,t)/3+soli.q(0,0,t).^2+soli.G),0,soli.t(Nti));
            end

            % function that calculates theta at time tstart for all (x,y)
%             soli.th.sa = sqrt(soli.a(soli.X,soli.Y,tstart)/12);
            soli.th.f = @(x,y) interp2(soli.x,soli.y,(soli.th.intx + soli.th.inty),x,y,'spline') + soli.th.intt(end);
            
     % Initial soliton approximation
        soli.us = @(x,y) soli.u(soli.th.f(x-soli.x0,y-soli.y0),x-soli.x0,y-soli.y0,tstart,soli.a,soli.G);
        
     % Assume can get with asymptotic approximation of 0
     soli.ua = @(x,y,t) zeros(size(x));
     soli.uax =@(x,y,t) zeros(size(x));
     soli.uay =@(x,y,t) zeros(size(x));

        end
