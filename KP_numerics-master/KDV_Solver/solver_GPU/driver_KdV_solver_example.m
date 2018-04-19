%% Driver file for KP2 numerical solver
%% Domain, ICs, etc go here
save_on  = 1;  % Set to nonzero if you want to run the solver, set
               % to 0 if you want to plot
rmdir_on = 0;  % Set to nonzero if you want to delete and remake the chosen directory
               % Useful for debugging
gpu_on   = 1;  % set to nonzero to use GPU, otherwise CPU
periodic = 1;  % set to nonzero to run periodic solver (no BCs need)
               % set to 0 to run solver with time-dependent BCs
plot_on  = 0;  % Set to 1 if you want to plot just before and just
               % after (possibly) calling the solver
check_IC = 0;  % Set to nonzero to plot the ICs and BCs without running the solver

for ii = 1
    %% Numerical Parameters
    tmax   = 10;      % Solver will run from t=0 to t = tmax
    numout = (tmax+1); % numout times will be saved (including ICs)
    Lx     = 300;     % Solver will run on x \in [-Lx,Lx]
    Nexp   = 9;
    Nx     = 2^Nexp;    % Number of Fourier modes in x-direction

    t      = linspace(0,tmax,numout);
   Nt      = 3;
   dt      = 10^(-Nt);
  %% Initial Condition and large-y approximation in time
  a = 2; x0 = 0;
  soli.ua = @(x) gpuArray(a*sech(sqrt(a/12)*(x-x0)).^2);
    %% Initial condition
  	soli.u0 = @(x)    soli.ua(x);
  
      ic_type = ['_soli_',...
                      '_a_',num2str(a),...
                      '_x0_',num2str(x0)];

    %% Generate directory, save parameters
	  q = strsplit(pwd,filesep);
    %% use CPU or GPU code, depending
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
    
    %% Choose directory data will be saved to, depending on computer
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
        
    %% Create directory run will be saved to
    data_dir = [maindir,slant,'Numerics',slant,'KP',slant,...
                'KdV',...
                '_tmax_',   num2str(round(tmax)),...
                '_Lx_',     num2str(Lx),...
                '_Nx_',     num2str(Nx),...
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
%     dd.(['qm',num2str(find(qrs==ql))]) = data_dir;
    savefile = sprintf('%sparameters.mat',data_dir);

    %% If chosen, run the solver using the parameters and conditions above
    if save_on
        if plot_on
            %% Load initial data
              xplot  = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
              tplot  = linspace(0,tmax,floor(tmax*10));
              u_init = soli.u0(xplot);
            %% Plot initial conditions and boundary conditions
            fontsize = 12;
            figure(1); clf;
                plot(xplot,u_init);
                title('Initial Conditions');
                xlabel('x'); ylabel('u');
            set(gca,'fontsize',fontsize,'fontname','times');
            pause(0.25);
            if check_IC
                drawnow;  pause(0.5);
                continue;
            end
        end
        % Save parameters
            save(savefile,'t','Nx','Lx',...
                              'soli','Nt',...
                              'periodic');
        try
        % Run KdV to determine boundary solver
            KdV_solver_periodic(  t, Lx, Nx, Nt, soli.u0, data_dir);
        catch e
            %send_mail_message('mdmaide2','Matlab',['Simulation failed with error ',e.message])
            disp(['Simulation failed with error ',e.message]);
	          continue;
        end
    else
        load(savefile);
    end


    if plot_on
        try
            plot_data_fun(data_dir);
            figure(4);
            print('sim','-dpng');
            send_mail_message('mdmaide2','Matlab',['Simulation ',data_dir,'done'],'sim.png')
        end
    else
        try
            send_mail_message('mdmaide2','Matlab',['Simulation ',data_dir,'done'])
        end
    end
end
