%% Driver file for KP2 numerical solver
%% Domain, ICs, etc go here
save_on  = 1;  % Set to nonzero if you want to run the solver, set
               % to 0 if you want to plot
periodic = 1;  % setto nonzero to run periodic solver (no BCs need)
               % set to 0 to run solver with time-dependent BCs                
plot_on  = 1;  % Set to 1 if you want to plot just before and just
               % after (possibly) calling the solver          
check_IC = 0;  % Set to nonzero to plot the ICs and BCs without running the solver

%% Numerical Parameters
tmax   = 5;      % Solver will run from t=0 to t = tmax
numout = 10*tmax+1; % numout times will be saved (including ICs)
Lx     = 50;     % Solver will run on x \in [-Lx,Lx]
Ly     = 30;     % Solver will run on y \in [-Ly,Ly]
Nx     = 2^8;    % Number of Fourier modes in x-direction
Ny     = 2^8;    % Number of Fourier modes in y-direction

t      = linspace(0,tmax,numout);

%% Initial Condition and large-y approximation in time
    ic_type = 'RP_soli';
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

    
   %% Two-Soliton amplitude, angle parameters
%     sam  = 2;
%     qm   = 1 ;
%     sap  = 1.5;
%     qp = qm - sam + sap;
% %     % Initial condition    
% %     [u0x,u0y,t0,u0mat] = RW1wave_evolve_fun(sam,qm,sap,qp,-Lx,Lx,-Ly,Ly,1);
% %     u0   = @(x,y) interp2(u0x,u0y,u0mat,x,y);
%     % Asymptotic approximation
%     aa     = @(x,t) sam^2 .* (x<0) + sap^2 .* (x>=0);
%     qa     = @(x,t)  qm   .* (x<0) +  qp   .* (x>=0);
%     uasoli = @(theta,x,t,a) a(x,t).*(sech(sqrt(a(x,t)/12).*theta)).^2;
%     theta = @(x,y,t,a,q) (x + q(x,t).*y - (a(x,t)/3+q(x,t).^2) * t);
%     uasy  = @(x,y,t) uasoli(theta(x,y,t,aa,qa),x,t,aa);
%     u0    = @(x,y)   uasy(x,y,0);
%% Generate directory, save parameters
if strcmp(computer,'MACI64')
    maindir = '/Volumes/Data Storage/Numerics/KP';
else
    maindir = 'H:';
end
    slant = filesep;
    
%% Create directory run will be saved to
data_dir = [maindir,slant,'Numerics',slant,'KP',slant,...
            '_tmax_',   num2str(round(tmax)),...
            '_Lx_',     num2str(Lx),...
            '_Nx_',     num2str(Nx),...
            '_Ly_',     num2str(Ly),...
            '_Ny_',     num2str(Ny),...
            '_bndry_condns_periodic',...
            '_init_condns_',ic_type,...
            slant];
% Create the data directory if necessary
if ~exist(data_dir,'dir')
    mkdir(data_dir);
else
    disp(['Warning, directory ',data_dir]);
    disp('already exists, possibly overwriting data');
end

savefile = sprintf('%sparameters.mat',data_dir);

%% If chosen, run the solver using the parameters and conditions above
if save_on
    % Load initial data
      xplot  = (2*Lx/Nx)*[-Nx/2:Nx/2];
      yplot  = (2*Ly/Ny)*[-Ny/2:Ny/2];
      [XPLOT,YPLOT] = meshgrid(xplot,yplot);
      tplot  = linspace(0,tmax,floor(tmax*10));
      u_init = u0(XPLOT,YPLOT);
    if plot_on
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
        pause(0.25);
        if check_IC
%             colorbar;
%             legend(ic_type);
            drawnow; 
            return;
        end
    end
    
    % Save parameters
        save(savefile,'t','Nx','Lx',...
                          'Ny','Ly',...
                          'u0','uasy','dxuasy','dyuasy',...
                          'periodic');
    % Run timestepper
        KP_solver_periodic( t, Lx, Nx,...
                               Ly, Ny,...
                               u0, uasy, dxuasy,dyuasy,...
                               data_dir );     
else
    load(savefile);
end

if plot_on
    plot_data_fun_2D(data_dir);
end
