%% Wrapper script for making experimental publication plots
% close all;
if strcmp(computer,'MACI64')
    updir = '/Users/dhl/Documents/MATLAB/KP_Whitham/';
    slant = '/';
else
    updir = 'D:\';
    slant = '\';
end

addpath(updir);

loadmat = 1; % Set to 1 to load in edge data, 0 to reformat data
runfigs = [1 0 1 0]; % [DSW_tunn, RW_tunn, RW_trap, DSW_trap]

if runfigs(1)
    % 1-wave
    figname = '1wave_evolve';
        xmin = -200;   % Set to 0 to begin as close to nozzle as possible
        xmax = 0; % Set to Inf to maximize upper conduit
        ymin = -300;   % Set to 0 to use first time
        ymax = 0; % Set to Inf to use last time 
        t    = 100;
        [x,y,t,U] = RW1wave_evolve_fun(xmin,xmax,ymin,ymax,t);
        cmin = 0;
        cred = 0;
        cmax = max(U(:)); % Set to Inf to maximize conduit area
            fhname = [pwd,slant,figname,'_contour'];
            [fh]   = make_contours_for_pub(U,x,y,ymin,ymax,...
                            xmin,xmax,cmin,cred,cmax,...
                            1,fhname,slant);
end

if runfigs(2)
    % 1-wave IC
    figname = '1wave_IC';
        xmin = -100;   % Set to 0 to begin as close to nozzle as possible
        xmax = 100; % Set to Inf to maximize upper conduit
        ymin = -75;   % Set to 0 to use first time
        ymax = 75; % Set to Inf to use last time 
        t    = 0.1;
        [x,y,t,U] = RW1wave_evolve_fun(xmin,xmax,ymin,ymax,t);
        cmin = 0;
        cred = 0;
        cmax = max(U(:)); % Set to Inf to maximize conduit area
            fhname = [pwd,slant,figname,'_contour'];
            [fh]   = make_contours_for_pub(U,x,y,ymin,ymax,...
                            xmin,xmax,cmin,cred,cmax,...
                            1,fhname,slant);
end


if runfigs(3)
    % 2-wave
    figname = '2wave_evolve';
        xmin = -50;   % Set to 0 to begin as close to nozzle as possible
        xmax =   75; % Set to Inf to maximize upper conduit
        ymin = - 100;   % Set to 0 to use first time
        ymax =   100; % Set to Inf to use last time 
        t    =   100;
        [x,y,t,U] = RW2wave_evolve_fun(xmin,xmax,ymin,ymax,t);
        cmin = 0;
        cred = 0;
        cmax = max(U(:)); % Set to Inf to maximize conduit area
            fhname = [pwd,slant,figname,'_contour'];
            [fh]   = make_contours_for_pub(U,x,y,ymin,ymax,...
                            xmin,xmax,cmin,cred,cmax,...
                            1,fhname,slant);
end


if runfigs(4)
    % 2-wave IC
    figname = '2wave_IC';
        xmin = -50;   % Set to 0 to begin as close to nozzle as possible
        xmax =  50; % Set to Inf to maximize upper conduit
        ymin = -100;   % Set to 0 to use first time
        ymax =  100; % Set to Inf to use last time 
        t    = 1;
        [x,y,t,U] = RW2wave_evolve_fun(xmin,xmax,ymin,ymax,t);
        cmin = 0;
        cred = 0;
        cmax = max(U(:)); % Set to Inf to maximize conduit area
            fhname = [pwd,slant,figname,'_contour'];
            [fh]   = make_contours_for_pub(U,x,y,ymin,ymax,...
                            xmin,xmax,cmin,cred,cmax,...
                            1,fhname,slant);
end