function KdV_solver_periodic(  t, Lx, Nx, Nt, u0, data_dir );
% Solves: KP eq. u_t + uu_x + u_xxx = 0
% on [-Lx,Lx] by FFT in space with
% Pseudospectral in space
% integrating factor v = exp[+i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat
% and RK4 in time
%
% Inputs:
%
% t        :  1D array of output times
% Lx, Nx   :  spatial grid parameters (x)
% u0       :  Initial condition, as a function
% data_dir : output directory where data at each output timestep
%             is written to a file #####.mat
%
% Outputs:  NONE except for data written to files

%% Set up temporal vector, increment counter, data directory
global tout inc dir
    tout = gpuArray([t(1):(10^(-Nt)):t(end)]);
    inc  = 0;
    dir  = data_dir;
    dt   = gpuArray(10^-Nt);

%% Define spatial and wavenumber grids
domain = struct;
    domain.x = gpuArray((2*Lx/Nx)*(-Nx/2:Nx/2-1)');
   domain.dx = gpuArray((2*Lx/Nx));
   domain.Lx = gpuArray(Lx);
   domain.kx = gpuArray((pi/Lx)*[0:Nx/2-1 0 -Nx/2+1:-1]');
    
%% Save variables, prepping to call the solver
    % Output what we are about to do
    disp(['Solving KdV eqtn.']);
    disp([' Time interval:  [', num2str(t(1)),...
           ',',num2str(t(end)),']']);
    % Construct initial condition on spatial domain
    u_init = gpuArray(u0(domain.x));
    u      = u_init; tnow = min(t);
    
    % Save initial condition (inc = 0)
	save(strcat(data_dir,num2str(inc,'%05d')),...
                'u','u_init',...
                'inc','tnow');
	inc = inc + 1;
    disp(['Time = ',num2str(tout(inc)),', inc = ',int2str(inc),...
          '/',int2str(length(tout))]);
	start = tic;

%% Call the solver
    myRK4_KdV( u_init, dt, domain)
 
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
   save([data_dir,'parameters.mat'],'finish','-append');
