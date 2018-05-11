%% Driver file for KP2 numerical solver
%% Domain, ICs, etc go here
rip;
save_on  = 1;  % Set to nonzero if you want to run the solver, set
               % to 0 if you want to plot
rmdir_on = 0;  % Set to nonzero if you want to delete and remake the chosen directory
               % Useful for debugging
gpu_on   = 1;  % set to nonzero to use GPU, otherwise CPU
periodic = 1;  % set to nonzero to run periodic solver (no BCs need)
               % set to 0 to run solver with time-dependent BCs
plot_on  = 1;  % Set to 1 if you want to plot just before and just
               % after (possibly) calling the solver
check_IC = 1;  % Set to nonzero to plot the ICs and BCs without running the solver

dd = struct();
sads   = 2;
qads   = 0;
x0s    = 100;
  %% Numerical Parameters
    tmax   = 10;      % Solver will run from t=0 to t = tmax
    numout = (tmax+1); % numout times will be saved (including ICs)
    Lx     = 300;     % Solver will run on x \in [-Lx,Lx]
    Ly     = Lx;     % Solver will run on y \in [-Ly,Ly]
    Nexp   = 8;
    Nx     = 2^Nexp;    % Number of Fourier modes in x-direction
    Ny     = 2^(Nexp);    % Number of Fourier modes in y-direction

    t      = linspace(0,tmax,numout);
   Nt      = 3;
   dt      = 10^(-Nt);
   
  %% Initial Condition and large-y approximation in time
   % One-soliton
    sau = 1; sad = sads;
    qau = 0; qad = qads;
    x0 = x0s; y0 = 0; x0odd = -x0s;
    w = Ly*2/3; tanw = 1/20; tstart = 30;
   % Define Soliton parameters and derivatives in x and y
      [ soli ] = vertical_segment(sau^2,sad^2,sad^2,...
                                     qau,qad,qad,...
                                     x0,y0,w,Lx,Ly,Nx,Ny,tstart);
      
      % Correction to zero mean; upper and lower bounds for KdV solver
        % Calculate dip amplitude for all y
        [soli] = zero_mean(soli,Ly,Lx,Ny,x0odd);   	
    	soli.hi = integral(@(x)soli.dip(x,-Ly,0),-Lx,+Lx);
    	soli.lo = integral(@(x)soli.dip(x,+Ly,0),-Lx,+Lx);
    	soli.ubhi  = @(x) -soli.hi/2*sech(x).^2;
    	soli.ublo  = @(x) -soli.lo/2*sech(x).^2;
    
      ic_type = ['_solikink_',...
                      '_sau_',num2str(sau),'_qu_',num2str(qau),...
                      '_sad_',num2str(sad),'_qd_',num2str(qad),...
                      '_x0_', num2str(x0 ),'_y0_',num2str(y0 )];

    %% Generate directory, save parameters
	  q = strsplit(pwd,filesep);
     % use CPU or GPU code, depending
        cpath = [strjoin(q(1:end-1),filesep),filesep,'solver'];
        gpath = [strjoin(q(1:end-1),filesep),filesep,'solver_GPU'];
        pathCell = regexp(path, pathsep, 'split');
    
        if gpu_on
        addpath(gpath)
        if any(strcmpi(cpath, regexp(path, pathsep, 'split')))
            rmpath(cpath)
        end
    else
        addpath(cpath)
        if any(strcmpi(gpath, regexp(path, pathsep, 'split')))
            rmpath(gpath)
        end
        end
    
     % Choose directory data will be saved to, depending on computer
        if strcmp(computer,'MACI64')
            maindir = '/Volumes/Data Storage/Numerics/KP';
            if ~exist(maindir,'dir')
                q = strsplit(pwd,filesep);
                maindir = strjoin(q(1:end-1),filesep);
            end
        elseif strcmp(computer,'PCWIN64')
            maindir = 'H:';
        else
            q = strsplit(pwd,filesep);
            maindir = '/rc_scratch/mima6446';%strjoin(q(1:end-1),filesep);
        end
    
        slant = filesep;
        
     % Create directory run will be saved to
        data_dir = [maindir,slant,'Numerics',slant,'KP',slant,...
                    '_tmax_',   num2str(round(tmax)),...
                    '_Lx_',     num2str(Lx),...
                    '_Nx_',     num2str(Nx),...
                    '_Ly_',     num2str(Ly),...
                    '_Ny_',     num2str(Ny),...
                    '_bndry_condns_periodic',...
                    '_init_condns_',ic_type,...
                    slant];
        dd.(['r',num2str(ii)])=data_dir;
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
%     dd.(['qm',num2str(find(qrs==ql))]) = data_dir;
    savefile = sprintf('%sparameters.mat',data_dir);

    %% If chosen, run the solver using the parameters and conditions above
    if save_on
        %% Generate zero-mean BCs from KdV
        [soli] = KdV_BC_generator( soli, t, Lx, Nx, Nt, maindir );

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
                cmap = load('CoolWarmFloat257.csv');
                [ cmap ] = cmap_edit( 0, gather(min(u_init(:))), gather(max(u_init(:))), cmap );
                contourf(XPLOT,YPLOT,u_init,length(cmap),'edgecolor','none'); xlabel('x'); ylabel('y');
                title('Initial Conditions');
                colormap(cmap);
                ax1 = gca;
            subplot(2,2,2)
                contourf(XPLOT,YPLOT,soli.ua(XPLOT,YPLOT,0),length(cmap),'edgecolor','none'); xlabel('x'); ylabel('y');
                title('Asymptotic u');
                set(gca,'CLim',get(ax1,'CLim'));
                colormap(cmap);
            subplot(2,2,3)
                contourf(XPLOT,YPLOT,soli.uax(XPLOT,YPLOT,0),100,'edgecolor','none'); xlabel('x'); ylabel('y');
                title('Asymptotic u, x-deriv');
                set(gca,'CLim',get(ax1,'CLim'));
                colormap(cmap);
            subplot(2,2,4)
                contourf(XPLOT,YPLOT,soli.uay(XPLOT,YPLOT,0),100,'edgecolor','none'); xlabel('x'); ylabel('y');
                title('Asymptotic u, y-deriv');
                set(gca,'CLim',get(ax1,'CLim'));
                colormap(cmap);

            set(gca,'fontsize',fontsize,'fontname','times');
            pause(0.25);
            if check_IC
    %             colorbar;
    %             legend(ic_type);
                drawnow;  pause(0.5);
                return;
            end
        end
        %% Save parameters
            save(savefile,'t','Nx','Lx',...
                              'Ny','Ly',...
                              'soli','Nt',...
                              'periodic');
%         try
        
        %% Run timestepper
            KP_solver_periodic( t, Lx, Nx, Nt,...
                                   Ly, Ny,...
                                   soli,...
                                   data_dir );
%         catch e
%             %send_mail_message('mdmaide2','Matlab',['Simulation failed with error ',e.message])
%             disp(['KP Simulation failed with error ',e.message]);
% 	          continue;
%         end
    else
        load(savefile);
    end


    if plot_on
        try
            plot_data_fun_2D(data_dir);
            figure(4);
            print('sim','-dpng');
            send_mail_message('mdmaide2','Matlab',['Simulation ',data_dir,'done'],'sim.png')
        end
    else
        try
            send_mail_message('mdmaide2','Matlab',['Simulation ',data_dir,'done'])
        end
    end

