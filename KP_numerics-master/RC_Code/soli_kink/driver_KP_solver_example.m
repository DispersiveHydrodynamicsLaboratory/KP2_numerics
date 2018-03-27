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

dd = struct();
sads   = [1 2 2  2  1];
qads   = [1 1 0 -1 -1];
x0s    = [100 100 100 -100 -100];
for ii = 3:5
    %% Numerical Parameters
    tmax   = 110;      % Solver will run from t=0 to t = tmax
    numout = (tmax+1); % numout times will be saved (including ICs)
    Lx     = 300;     % Solver will run on x \in [-Lx,Lx]
    Ly     = Lx;     % Solver will run on y \in [-Ly,Ly]
    Nexp   = 9;
    Nx     = 2^Nexp;    % Number of Fourier modes in x-direction
    Ny     = 2^(Nexp);    % Number of Fourier modes in y-direction

    t      = linspace(0,tmax,numout);
   Nt      = 3;
   dt      = 10^(-Nt);
    %% Initial Condition and large-y approximation in time
%         ic_type = ['KP2_validation_Nexp_',num2str(Nexp),'_dt_',num2str(Nt),'_twosoli'];
        %% One-soliton
        sau = 1; sad = sads(ii);
        qau = 0; qad = qads(ii);
            x0 = x0s(ii); y0 = 0; x0odd = -x0s(ii); w = Ly*2/3;
            tanw = 1/20;
%% Define Soliton parameters and derivatives in x and y
    y     = (2*Ly/Ny)*(-Ny/2:Ny/2-1)';
    x     = (2*Lx/Nx)*(-Nx/2:Nx/2-1)';
    [X,Y] = meshgrid(x,y);
    au         = sau^2;
    ad         = sad^2;
    soli.Y     = gpuArray(Y);
    soli.tw    = gpuArray(tanw);
    soli.w     = gpuArray(w);
	soli.ybox  = 1/2.*(-tanh(soli.tw*(soli.Y-soli.w/2))+tanh(soli.tw*(soli.Y+soli.w/2)));
    soli.a     = @(x,y,t) gpuArray(( (ad-au)./2.*tanh(soli.tw.*(y)) + (ad+au)/2 )) ;
    soli.ax    = @(x,y,t) gpuArray(zeros(size(x)));
    soli.ay    = @(x,y,t) gpuArray((ad-au)/2*soli.tw*sech(soli.tw*(y)).^2);
    soli.q     = @(x,y,t) gpuArray((qad-qau)/2*tanh(soli.tw*(y)) + (qad+qau)/2);
    soli.qx    = @(x,y,t) gpuArray(zeros(size(x)));
    soli.qy    = @(x,y,t) gpuArray((qad-qau)/2*soli.tw*sech(soli.tw*(y)).^2);
    soli.x0    = gpuArray(x0);
    soli.y0    = gpuArray(y0);
    soli.G     = gpuArray(0);
%% Define Soliton function, derivatives
    soli.u  = @(theta,x,y,t,a,g) g + a(x,y,t).*(sech(sqrt(a(x,y,t)/12).*theta)).^2; % CORRECT
    soli.dx = @(theta,dxtheta,x,y,t,a,ax) sech(sqrt(a(x,y,t)/12).*theta).^2.*...
               (ax(x,y,t) + sqrt(a(x,y,t)/12).*tanh(sqrt(a(x,y,t)/12).*theta).*...
                    dxtheta);
    soli.dy = @(theta,dytheta,x,y,t,a,q,ay) sech(sqrt(a(x,y,t)/12).*theta).^2.*...
               (ay(x,y,t) + sqrt(a(x,y,t)/12).*tanh(sqrt(a(x,y,t)/12).*theta).*...
                    dytheta);
% Theta and derivatives
    soli.th = @(x,y,t,a,q,g) (x + q(x,y,t).*y - (a(x,y,t)/3+q(x,y,t).^2+g) * t); % FIXED
    soli.thx = @(x,y,t,a,ax,q,qx,g) ((g*t-x-y.*q(x,y,t)+t*q(x,y,t).^2).*ax(x,y,t)+...
                        a(x,y,t).*(-2+t*ax(x,y,t)-2*(y-2*t*q(x,y,t)).*qx(x,y,t)));
    soli.thy = @(x,y,t,a,ay,q,qy,g) ((g*t-x-y.*q(x,y,t)+t*q(x,y,t).^2).*ay(x,y,t)+...
                        a(x,y,t).*(t.*ay(x,y,t)-2*y.*qy(x,y,t)+...
                        q(x,y,t).*(-2+ 4*t.*qy(x,y,t))));
  
% Asymptotic soliton approximation
    soli.ua  = @(x,y,t)  soli.u(soli.th(x-soli.x0,y-soli.y0,t,soli.a,soli.q,soli.G),...
                                x-soli.x0,y-soli.y0,t,soli.a,soli.G);
    soli.uax  = @(x,y,t) soli.dx(soli.th(x-soli.x0,y-soli.y0,t,soli.a,soli.q,soli.G),...
                            soli.thx(x-soli.x0,y-soli.y0,t,soli.a,soli.ax,soli.q,soli.qx,soli.G),...
                             x-soli.x0,y-soli.y0,t,soli.a,soli.ax);
    soli.uay  = @(x,y,t) soli.dy(soli.th(x-soli.x0,y-soli.y0,t,soli.a,soli.q,soli.G),...
                                 soli.thy(x-soli.x0,y-soli.y0,t,soli.a,soli.ay,soli.q,soli.qy,soli.G),...
                                  x-soli.x0,y-soli.y0,t,soli.a,soli.q,soli.ay);
%% Initial condition, with odd reflection
	soli.x0odd = x0odd;
	soli.u0 = @(x,y)    soli.ua(x,y,0) - soli.ua(x+soli.x0-soli.x0odd,y,0).*soli.ybox;

    ic_type = ['_solikink_',...
                    '_sau_',num2str(sau),'_qu_',num2str(qau),...
                    '_sad_',num2str(sad),'_qd_',num2str(qad),...
                    '_x0_',num2str(x0 ),'_y0_',num2str(y0)];

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
        if plot_on
            %% Load initial data
              xplot  = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
              yplot  = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
              [XPLOT,YPLOT] = meshgrid(xplot,yplot);
              tplot  = linspace(0,tmax,floor(tmax*10));
              u_init = soli.u0(XPLOT,YPLOT);

            %% Plot initial conditions and boundary conditions
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
                continue;
            end
        end

        % Save parameters
            save(savefile,'t','Nx','Lx',...
                              'Ny','Ly',...
                              'soli','Nt',...
                              'periodic');
        try
        % Run timestepper
            KP_solver_periodic( t, Lx, Nx, Nt,...
                                   Ly, Ny,...
                                   soli,...
                                   data_dir );    
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
end
