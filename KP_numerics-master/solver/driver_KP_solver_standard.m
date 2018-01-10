% Driver script for running the KP equation solver
% KP_solver_periodic.m
save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 1; % setget to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs                
plot_on  = 0;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver          
check_IC = 0; % Set to nonzero to plot the ICs and BCs without running the solver

%% Numerical Parameters
tmax     = 5;    % Solver will run from t=0 to t=tmax
xmax     = 100;   % Solver will solve on domain x=0 to z=xmax
                 % With an odd reflection of IC from x = -xmax to x=0
ymax     = 100;   % Solver will solve on domain y=-ymax to y=ymax
numout   = tmax+1;           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
% dxinit =  1/10;     % Spatial Discretization 
% dyinit =  1;     % KP Solver not valdiated yet                  
% Nx       = round(2*xmax/dxinit);
% Ny       = round(2*ymax/dyinit);
Nx = 2^9;
Ny = 2^9;
if periodic
    dx       = 2*xmax/Nx;    % Spatial  discretization
	dy       = 2*ymax/Ny;    % Spatial  discretization
else
    dx       = 2*xmax/(Nx+1);    % Spatial  discretization
    dy       = 2*ymax/(Ny+1);    % Spatial  discretization
end
%% KP PARAMETERS
epsilon = 1; % u_xxx scaling factor
lambda  = 1; % Sign on the mixed term (+1 is KP2, -1 is KP1)

%% PDE Initial and Boundary Conditions
%% Set up 'True' IC
%% Set up soliton function
a =@(x,y) 3; q =@(x,y) 1/4;
u = @(theta,x,y,a) a(x,y).*(sech(sqrt(a(x,y)/12).*theta)).^2;
theta = @(x,y,a,q,x0,y0) ( (x-x0) + q(x,y).*(y-y0) );
IC = @(X,Y) u(theta(X,Y,a,q,xmax/2,0),X,Y,a);

%% Factor zeroing IC as y->ymax
winspan = 2.5; % width of decay window
w = @(Y) repmat(tukeywin(size(Y,1),winspan/ymax),[1,size(Y,2)]);

%% Putting it all together
u0 = @(X,Y) (IC( X,Y)-IC(X+xmax,Y)).*w(Y);
    ic_type = '';
if strcmp(computer,'MACI64')
    maindir = '/Users/dhl/Documents/MATLAB';
    slant = '/';
else
    maindir = pwd;%'H:';
    slant = '\';
end

%% Create directory run will be saved to
data_dir = [maindir,slant,'Numerics',slant,'KP',slant,...
            '_tmax_',   num2str(round(tmax)),...
            '_xmax_',   num2str(round(xmax)),...
            '_Nx_',     num2str(Nx),...
            '_ymax_',   num2str(round(ymax)),...
            '_Ny_',     num2str(Ny),...
            '_epsilon_',num2str(epsilon),...
            '_lambda_', num2str(lambda),...
            '_bndry_condns_','periodic',...
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
      xplot  = dx*[-Nx/2:Nx/2];
      yplot  = dy*[-Ny/2:Ny/2];
      [XPLOT,YPLOT] = meshgrid(xplot,yplot);
      tplot  = linspace(0,tmax,floor(tmax*10));
      u_init = u0(XPLOT,YPLOT);
    if plot_on
        % Plot initial conditions and boundary conditions
        fontsize = 12;
        figure(1); clf;
        if ~periodic
            subplot(3,1,1);
        end
            contourf(XPLOT,YPLOT,u_init,100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
            title('Initial Conditions');
        set(gca,'fontsize',fontsize,'fontname','times');
        pause(0.25);
        if check_IC
            colorbar;
%             legend(ic_type);
            drawnow; 
            return;
        end
    end
    
    % Save parameters
        save(savefile,'t','Nx','dx','xmax',...
                          'Ny','dy','ymax',...
                          'u0','periodic','epsilon','lambda');
    % Run timestepper
        KP_solver_periodic( t, xmax, Nx,...
                               ymax, Ny,...
                               epsilon, lambda,...
                               u0, data_dir );     
else
    load(savefile);
end

if plot_on
    plot_data_fun_2D(data_dir);
end

% %% If chosen, plot data associated with the parameters and conditions above
% if plot_on
%     disp('Calculating maximum time increment in saved data files...');
%     for ii=1:length(t)+1
%         [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d.mat')),'r');
%         if fid == -1 % File does not exist
%             tm = ii-1;
%             disp(['Maximum time = ',num2str(t(tm))]);
%             break;
%         end
%         fclose(fid);
%     end
%     if ii == length(t)
%         disp(['Maximum time = ',num2str(t(tm))]);
%     end
%     % Get rid of larger t values
%     t = t(1:tm);
%     if periodic
%         xplot  = dx:dx:xmax;
%     else
%         xplot  = dx:dx:xmax-dx;
%     end
%     A_full = zeros(length(t)-1,length(xplot));
%     % Load first time step
%     load(strcat(data_dir,num2str(0,'%05d')),'A_init');
%         fontsize = 12;
%         fig=figure(2); clf;
%         plot(xplot,A_init);
%         qaxis = axis;
%         title(['Time: ',num2str(t(1))]);
%         set(gca,'fontsize',fontsize,'fontname','times');
%         drawnow
% %         input('Return');
%     %Plot subsequent time steps
%     for tind=2:tm
%         load(strcat(data_dir,num2str(tind,'%05d')),'A','tnow');
%         fontsize = 12;
%         fig=figure(2); clf;
%         plot(xplot,A);
%         A_full(tind-1,:) = A;
%         axis(qaxis)
%         hold off;
%         title(['Time: ',num2str(tnow)]);
%         set(gca,'fontsize',fontsize,'fontname','times');
%         drawnow;
%     end
% 
% end
