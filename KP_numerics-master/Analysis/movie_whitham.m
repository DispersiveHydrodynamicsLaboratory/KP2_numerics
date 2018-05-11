%% Script for generating movies out of various things
clear; close all; 
if strcmp(getenv('computername'),'CAKELIE2012')
    main_dir = 'H:\Numerics\KP\';
else
    main_dir = pwd;
end
run_dir  = '_tmax_150_Lx_800_Nx_1024_Ly_400_Ny_512_bndry_condns_periodic_init_condns__solitest_true_soln\';
data_dir = [main_dir, run_dir];
movie_name = 'true_soli_tstart_30_full';
% Parameters
pix = [1024]; % Number of pixels in horizontal direction
tfac = 12;  % Speedup in time
num_files = Inf; % set to number of pictures desired; if Inf, will use all pictures
fontsize = 20; % Fontsize for length and time scales
show_time = 1; % Set nonzero if you want display of current time

xlims = [-800 800];
ylims = [-400 400];
clims = [-0.1 1];

output_dir = [data_dir,'movie_files/'];

if ~exist(output_dir,'dir');
    mkdir(output_dir);
end

% Finds and saves axes parameters
load([data_dir,'parameters.mat'],'Nx','Ny','Lx','Ly','t');
%% x and y vectors; Check initial condition
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [X,Y] = meshgrid(x,y);
    [fh,axh]=plot_nice_contour(data_dir,0,1);
        set(axh,'XLim',xlims,'YLim',ylims,'CLim',clims,'fontsize',fontsize);
        set(fh,'Color','w'); 
        title(['T = 0']);
	disp(['Generating quicktime movie for horizontal pixel resolution ',...
            int2str(pix)]);
    % Create movie object
    movObj = QTWriter([output_dir,movie_name,'.mov'],...
                      'MovieFormat','Photo TIFF');%,'Quality',85);
	writeMovie(movObj,getframe(fh));

for ti=1:length(t)-1
    [fh,axh] = plot_nice_contour(data_dir,ti,1);
        set(axh,'XLim',xlims,'YLim',ylims,'CLim',clims,'fontsize',fontsize);
        set(fh,'Color','w'); 
    title(['T = ',num2str(ti)]);
        writeMovie(movObj,getframe(fh));
end
    % Default for movies to loop
    movObj.Loop = 'loop';

    disp(['Mean framerate = ',num2str(movObj.MeanFrameRate),' fps']);
    
    % Finish writing
    close(movObj);
    
    disp(['Movie can be found at ',output_dir,movie_name,'.mov']);
