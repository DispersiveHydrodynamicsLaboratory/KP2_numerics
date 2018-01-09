function KP_solver_periodic( t, xmax, Nx,...
                                ymax, Ny,...
                                epsilon, lambda,...
                                u0, directory )
% Solves: KP eq. (u_t + 6uu_x + epsilon^2 u_xxx)_x + lambda u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with integrating factor 
% v = exp[+i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat and RK4 in time
%
% Inputs:
%
% t        :  1D array of output times
% xmax, dx :  spatial grid parameters (x)
% ymax, dy :  spatial grid parameters (y)
% epsilon  :  scaling on u_xxx term
% lambda   :  scaling on mixed term
% u0       :  initial condition in space at t=0; function of X and Y
% directory - output directory where data at each output timestep
%             is written to a file #####.mat
%
% Outputs:  NONE except for data written to files
%
  global inc dir tout start;

  % Set global variables for ode solver output function
  tout = t;
  dir = directory;
  inc = 0;
  dt = 1e-3;

  % Setup grid
  Lx = xmax; Ly = ymax;
    x = (2*Lx/Nx)*(-Nx/2:Nx/2-1)';
    y = (2*Ly/Ny)*(-Ny/2:Ny/2-1)';
    [X,Y] = meshgrid(x,y);
    k = (pi/Lx)*[0:Nx/2-1 0 -Nx/2+1:-1]';
    l = (pi/Ly)*[0:Ny/2-1 0 -Ny/2+1:-1]';
    [KX,KY] = meshgrid(k,l);
    o = eps; % Attempt to remove possible issue of dividing by 0
    iphi = -1i*(-epsilon^2*KX.^3+lambda*KY.^2./(KX+1i*lambda*o));

%% Factor zeroing as y->ymax
winspan = 2.5; % width of decay window
w = @(Y) repmat(tukeywin(size(Y,1),winspan/Ly),[1,size(Y,2)]);
W = w(Y);
    
  %
  % Output what we are about to do
  %
  disp(['Solving KP eqtn.']);
  disp([' Time interval:  [', num2str(t(1)),...
        ',',num2str(t(end)),']']);
  
  %
  % Construct initial condition on spatial domain
  %
  u_init = u0(X,Y);
  u      = u_init; tnow = min(t);
  %
  % Construct initial condition on wavenumber domain
  %
  v_init = fft2(u0(X,Y));
  v      = v_init;
  % Save initial data 
  save(strcat(dir,num2str(inc,'%05d')),'u_init','v_init','inc');
  inc = inc + 1;
  save(strcat(dir,num2str(inc,'%05d')),'u','v','tnow','inc');
  disp(['Time = ',num2str(tout(inc)),', inc = ',int2str(inc),...
          '/',int2str(length(tout))]);
  
%   %Set appropriate ODE solver options.
%   opts = odeset('OutputFcn',@output,'Refine',1,'RelTol',1e-6, ...
%                 'AbsTol',1e-6,'Stats','on','maxstep',dt);
  start = tic;

        myRK4_KP2( iphi, KX,KY, W, dt, tout, v_init, dir)
 
  % Finish and clean up
  finish = toc(start);
  disp('Calculation Finished');
  days_left = datenum([0 0 0 0 0 toc(start)]);
  time_left=datevec(days_left-floor(days_left));
  disp(['Computation time = ',...
        int2str(floor(days_left)),'d ',...
        int2str(time_left(4)),'h ',...
        int2str(time_left(5)),'m ',...
        num2str(time_left(6)),'s']);
  
%
% % Called by odesolver at every output time
% %
% function status = output(t,v,flag)
%   global inc dir tout start;
%   status = 0;
%   if strcmp(flag,'')
%     for ii=1:length(t)
%       % Extract solution
%       v = v(:,ii);
%       tnow = t(ii);
%       %%%%% TODO: CONVERT u TO ifft2[exp[-iphi*t]*v]
%       u = nans(size(v)); 
%       % Increment and output estimate of time left
%       inc = inc + 1;
%       % Output results to file for subsequent analysis
%       save(strcat(dir,num2str(inc,'%05d')),'u','v','tnow','inc');
%       percentage_of_work = (inc)/(length(tout));
%       seconds_left = ((1/percentage_of_work)-1)*(toc(start));
%       days_left = datenum([0 0 0 0 0 seconds_left]);
%       time_left=datevec(days_left-floor(days_left));
%       disp(['Time = ',num2str(tout(inc)),...
%           ', inc = ',int2str(inc),'/',int2str(length(tout)),...
%           ', computation time left = ',...
%           int2str(floor(days_left)),'d ',...
%           int2str(time_left(4)),'h ',...
%           int2str(time_left(5)),'m ',...
%           num2str(time_left(6)),'s']);
%     end
%   end
  




  