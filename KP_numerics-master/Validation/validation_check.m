%% Processes validation files
this_dir = regexp(pwd,filesep,'split');
% sol_dir  = [this_dir(1:end-1),
plot_on = 1;
save_on = 0;

Nexps = 6:8;
Ne_max_err = zeros(1,length(Nexps));

if save_on
    uNx = struct();
end

% for tnow = 10;% [9 10 11]
    for Ni = 1:length(Nexps)
        %% Numerical Parameters
        tmax   = 5;      % Solver will run from t=0 to t = tmax
        numout = 10*tmax+1; % numout times will be saved (including ICs)
        Lx     = 50;     % Solver will run on x \in [-Lx,Lx]
        Ly     = Lx/2;     % Solver will run on y \in [-Ly,Ly]
        Nx     = 2^Nexps(Ni);    % Number of Fourier modes in x-direction
        Ny     = 2^Nexps(Ni)/2;    % Number of Fourier modes in y-direction

    t      = linspace(0,tmax,numout);
   Nt      = 4;
   dt      = 10^(-Nt);

        %% Initial Condition and large-y approximation in time
%             ic_type = ['KP2_validation_Nexp_',num2str(Nexps(Ni))];
%         ic_type = ['KP2_validation_Nexp_',num2str(Nexp),'_dt_',num2str(4)];
        ic_type = ['KP2_validation_Nexp_',num2str(Nexps(Ni)),'_dt_',num2str(Nt),'_straightsoli'];

            %% One-soliton, corrected for nonzero integral in x
                sa = 0.5; q = 0; x0 = 0;
    %             [ soli ] = one_soli(sa,q,x0,Lx);

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
                disp(data_dir);
        %% Check if dir exists
        if exist(data_dir,'dir')
            disp(['Nexp ',num2str(Nexps(Ni)),' exists']);
        else
            disp(['Nexp ',num2str(Nexps(Ni)),' is missing']);
            continue;
        end
        %% Load final time, evolved IC
        load([data_dir,'parameters.mat'],'soli');
        load([data_dir,num2str(25,'%05d'),'.mat'],'u','tnow');
        %% Compare 
        x = (2*Lx/Nx)*(-Nx/2:Nx/2-1)';
        y = (2*Ly/Ny)*(-Ny/2:Ny/2-1)';
        [X,Y] = meshgrid(x,y);
        UA = soli.ua(X,Y,tnow);
        Ne_max_err(Ni) = max(max(abs(UA-u)));
        if save_on
            uNx.(['u',num2str(Nexps(Ni))]) = u;
            uNx.(['UA',num2str(Nexps(Ni))]) = UA;
        end

        %% Comparison figure
        if plot_on
            figure(2); clf;
                contourf(X,Y,UA,100,'edgecolor','none');
                colorbar; title(['tnow: ',num2str(tnow)]);
            figure(3); clf;
                contourf(X,Y,u,100,'edgecolor','none');
                colorbar; title(['tnow: 10']);
                drawnow; % input('R');
        end

        disp('');
    end
% end
figure(1); clf;
loglog(2.^(Nexps),(Ne_max_err),'x');

if save_on
	save('Validation_data.mat','uNx','Nexps','Ne_max_err');
end
