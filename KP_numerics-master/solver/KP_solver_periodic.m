function KP_solver_periodic( t, Lx, Nx,...
                                Ly, Ny,...
                                u0, uasy,...
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

%% Set up temporal vector, increment counter
global tout inc
    tout = t;
    inc  = 0;
    dt   = 10^-3;

%% Rescale spatial grid to be from [-pi,pi]
    x = (2*pi/Nx)*(-Nx/2:Nx/2-1)';
    y = (2*pi/Ny)*(-Ny/2:Ny/2-1)';
    [X,Y] = meshgrid(x,y);
    k = [0:Nx/2-1 0 -Nx/2+1:-1]';
    l = [0:Ny/2-1 0 -Ny/2+1:-1]';
    [KX,KY] = meshgrid(k,l);
    u0 = @(x,y) u0(Lx/pi*x,Ly/pi*y);
    uasy = @(x,y,t) uasy(Lx/pi*x,Ly/pi*y,t);
    
%% Windowing function (from Kao 2010), rescaled, and derivs
n = 27; an = (1.111)^n*log(10);
W = @(y) exp( an * abs(y/pi)^n );
Wp = @(y) an*n/p * W(y) * abs(y/pi)^(n-1)*sign(y/pi);
Wpp = @(y) an*n/p^2 * W(y) * abs(y/pi)^(n-2) * ...
           ( (n-1 + an*n*abs(y/pi)^n)*sgn(y/pi)^2 );
    
%% Write static matrices here, so not constantly redefining
    o = eps; % Attempt to remove possible issue of dividing by 0
    % integrating factor exponent
    iphi = 1i*(pi*Lax/(ly^2)*KY.^2./(KX+1i*o)-(pi/Lx)^3*KX.^3); 
    %% MM: NEEDS WORK: NEEDS D1x and D1y for 2D 
    % First derivative (2nd order centered diff)
	D1  = 1/(2*dz) * spdiags([ -ones(Nz,1),...
                                  zeros(Nz,1),...
                                  ones(Nz,1)],...
                               -1:1, Nz,Nz);
	% Second Derivative (2nd order centered diff)
	D2  = 1/(dz^2) * spdiags([  ones(Nz,1),...
                               -2*ones(Nz,1),...
                                  ones(Nz,1) ],...
                               -1:1,Nz,Nz);

%% Save variables, prepping to call the solver
    % Output what we are about to do
    disp(['Solving KP eqtn.']);
    disp([' Time interval:  [', num2str(t(1)),...
           ',',num2str(t(end)),']']);
    % Construct initial condition on spatial domain
    u_init = u0(X,Y);
    u      = u_init; tnow = min(t);
    % Construct initial v, the windowed version of u
    v_init = u_init - (1 - W(Y)) * uasy(X,Y,0);
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
    myRK4_KP2( Vhat_init, uasy, dt, tout, W, Wp, Wpp,...
                    iphi, D1x, D1y, KX)
 
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
    
