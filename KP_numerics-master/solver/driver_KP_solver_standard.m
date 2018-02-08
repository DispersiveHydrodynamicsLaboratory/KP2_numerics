%% Driver file for KP2 numerical solver
%% Domain, ICs, etc go here
save_on  = 1;  % Set to nonzero if you want to run the solver, set
               % to 0 if you want to plot
rmdir_on = 1;  % Set to nonzero if you want to delete and remake the chosen directory
               % Useful for debugging
periodic = 1;  % setto nonzero to run periodic solver (no BCs need)
               % set to 0 to run solver with time-dependent BCs                
plot_on  = 0;  % Set to 1 if you want to plot just before and just
               % after (possibly) calling the solver          
check_IC = 0;  % Set to nonzero to plot the ICs and BCs without running the solver

    %% Numerical Parameters
    tmax   = 10;      % Solver will run from t=0 to t = tmax
    numout = 5*tmax+1; % numout times will be saved (including ICs)
    Lx     = 50;     % Solver will run on x \in [-Lx,Lx]
    Ly     = Lx/2;     % Solver will run on y \in [-Ly,Ly]
    Nx     = 2^9;    % Number of Fourier modes in x-direction
    Ny     = Nx/2;    % Number of Fourier modes in y-direction

    t      = linspace(0,tmax,numout);

    %% Initial Condition and large-y approximation in time
        ic_type = ['KP2_example_',num2str(Nexp)];
        %% One-soliton, corrected for nonzero integral in x
            sa = 0.5; q = 0.75; x0 = 0;
            [ soli ] = one_soli(sa,q,x0,Lx);

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
        if rmdir_on == 1
            disp('already exists, rewriting entire folder');
            rmdir(data_dir,'s')
            mkdir(data_dir);
        else
            disp('already exists, possibly overwriting data');
        end
    end

    savefile = sprintf('%sparameters.mat',data_dir);

    %% If chosen, run the solver using the parameters and conditions above
    if save_on
        if plot_on
            % Load initial data
              xplot  = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
              yplot  = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
              [XPLOT,YPLOT] = meshgrid(xplot,yplot);
              tplot  = linspace(0,tmax,floor(tmax*10));
              u_init = soli.u0(XPLOT,YPLOT);

            % Plot initial conditions and boundary conditions
            fontsize = 12;
            figure(1); clf;
            subplot(2,2,1)
                contourf(XPLOT,YPLOT,u_init,100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Initial Conditions');
            subplot(2,2,2)
                contourf(XPLOT,YPLOT,soli.ua(XPLOT,YPLOT,0),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Asymptotic u');
            subplot(2,2,3)
                contourf(XPLOT,YPLOT,soli.uax(XPLOT,YPLOT,0),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Asymptotic u, x-deriv');
            subplot(2,2,4)
                contourf(XPLOT,YPLOT,soli.uay(XPLOT,YPLOT,0),100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
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
                              'soli',...
                              'periodic');
        % Run timestepper
            KP_solver_periodic( t, Lx, Nx,...
                                   Ly, Ny,...
                                   soli,...
                                   data_dir );     
    else
        load(savefile);
    end
end
if plot_on
    plot_data_fun_2D(data_dir);
end
