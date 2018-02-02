function KP_solver_periodic( t, Lx, Nx,...
                                Ly, Ny,...
                                u0, uasy, dxuasy, dyuasy,...
                                data_dir )
% Solves: KP eq. (u_t + uu_x + u_xxx)_x + u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with
% Window technique described in Kao, Kodama 2010 in y-direction
% Pseudospectral and 2nd-order centered differences in space
% integrating factor v = exp[+i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat 
% and RK4 in time
%
% Inputs:
%
% t        :  1D array of output times
% Lx, Nx   :  spatial grid parameters (x)
% Ly, Ny   :  spatial grid parameters (y)
% u0       :  initial condition in space at t=0; function of X and Y
% uasy     :  what u would look like in the far field without the window
% data_dir : output directory where data at each output timestep
%             is written to a file #####.mat
%
% Outputs:  NONE except for data written to files

%% Set up temporal vector, increment counter, data directory
global tout inc dir
    tout = t;
    inc  = 0;
    dir  = data_dir;
    dt   = 10^-3;

%% Rescale spatial grid to be from [-pi,pi]
domain = struct;
    domain.x = (2*pi/Nx)*(-Nx/2:Nx/2-1)';
   domain.dx = (2*pi/Nx);
   domain.Lx = Lx;
    domain.y = (2*pi/Ny)*(-Ny/2:Ny/2-1)';
   domain.dy = (2*pi/Ny);
   domain.Ly = Ly;
    [domain.X,domain.Y] = meshgrid(domain.x,domain.y);
    domain.k = [0:Nx/2-1 0 -Nx/2+1:-1]';
    domain.l = [0:Ny/2-1 0 -Ny/2+1:-1]';
    [domain.KX,domain.KY] = meshgrid(domain.k,domain.l);
    % Rescale IC, asymptotic solution accordingly
        u0   = @(x,y) u0(Lx/pi*x,Ly/pi*y);
        uasy = @(x,y,t)   uasy(Lx/pi*x,Ly/pi*y,t);
      dxuasy = @(x,y,t) dxuasy(Lx/pi*x,Ly/pi*y,t);
      dyuasy = @(x,y,t) dyuasy(Lx/pi*x,Ly/pi*y,t);
    
%% Write static matrices here, so not constantly redefining  
    %% Windowing function (from Kao 2010), rescaled, and derivs
    n = 27; an = (1.111)^n*log(10);
    W = exp( -an * abs(domain.Y/pi).^n );
    Wp = -an*n/pi * W .* abs(domain.Y/pi).^(n-1).*sign(domain.Y/pi);
    Wpp = an*n/pi^2 * W .* abs(domain.Y/pi).^(n-2) .* ...
               ( (-(n-1) + an*n*abs(domain.Y/pi).^n).*sign(domain.Y/pi).^2 );
	%% iphi
    o = eps; % Attempt to remove possible issue of dividing by 0
    % integrating factor exponent
    iphi = 1i*(pi*Lx/(Ly^2)*domain.KY.^2./(domain.KX+o)-(pi/Lx)^3*domain.KX.^3); 
    
%% Save variables, prepping to call the solver
    % Output what we are about to do
    disp(['Solving KP eqtn.']);
    disp([' Time interval:  [', num2str(t(1)),...
           ',',num2str(t(end)),']']);
    % Construct initial condition on spatial domain
    u_init = u0(domain.X,domain.Y);
    u      = u_init; tnow = min(t);
    % Construct initial v, the windowed version of u
    v_init = u_init - (1 - W) .* uasy(domain.X,domain.Y,0);
    v      = v_init;
    % Construct initial Vhat, the fft'd, integrating factored v
    Vhat_init = fft2(v_init);
    Vhat      = Vhat_init;
    
    % Save initial condition (inc = 0)
	save(strcat(data_dir,num2str(inc,'%05d')),...
                'u','u_init','v','v_init','Vhat',...
                'Vhat_init','inc','tnow');
	inc = inc + 1;
    disp(['Time = ',num2str(tout(inc)),', inc = ',int2str(inc),...
          '/',int2str(length(tout))]);
	start = tic;

%% Call the solver
    myRK4_KP2( Vhat_init, uasy, dxuasy, dyuasy, dt, W, Wp, Wpp,...
                    iphi, domain )
 
%% Finish and clean up
    finish = toc(start);
    disp('Calculation Finished');
    days_left = datenum([0 0 0 0 0 toc(start)]);
    time_left=datevec(days_left-floor(days_left));
    disp(['Computation time = ',...
        int2str(floor(days_left)),'d ',...
        int2str(time_left(4)),'h ',...
        int2str(time_left(5)),'m ',...
        num2str(time_left(6)),'s']);
    
