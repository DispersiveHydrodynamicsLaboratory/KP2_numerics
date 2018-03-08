function[fig,ch] = plot_quick_contour(data_dir,varargin)
loaddir = data_dir;
q = strfind(loaddir,'camassa');
if q
    load([loaddir,'parameters.mat'],'t','Nx','dx','xmax'); %,'f');
    Nz = Nx; dz = dx; zmax = xmax;
    param = 'u';
else
    load([loaddir,'parameters.mat'],'t','Nz','dz','zmax'); %,'f');
    param = 'A';
end
z = dz*[1:Nz];
if nargin>1
    reformat_on = varargin{1};
    if nargin>2
        fignum = varargin{2};
        if nargin>3
            tmax = varargin{3};
        else
            tmax = Inf;
        end
    else
        fignum = 5;
        tmax = Inf;
    end
else
    fignum = 5;
    if exist([data_dir, 'conduit_edges.mat'],'file')
        reformat_on = 0;
    else
        reformat_on = 1;
    end
    tmax = Inf;
end
    parafile = sprintf('%sparameters.mat',data_dir);
    load(parafile);
    if reformat_on
        numfile = reformat_numerics_files(data_dir);
    else
        numfile = [data_dir, 'conduit_edges.mat'];
    end
        load(numfile);
    %% Set up grid to interpolate FROM
    if periodic
        zold  = dz:dz:zmax;
    else
        zold  = dz:dz:zmax-dz;
    end
    % Calculating maximum time increment in saved data files...
    for ii=1:length(t)+1
        [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d.mat')),'r');
        if fid == -1 % File does not exist
            tm = ii-1;
            break;
        end
        fclose(fid);
    end
    % Get rid of larger t values
        t = t(1:tm);
    %% Set up grid to interpolate TO
    N = 300;
    znew = linspace(0, zmax, N);
    tnew = linspace(min(t),min(tmax,t(tm)), N);
    
    area_new = interp2(zold,t,area_mat,znew,tnew');
    
    fig=figure(fignum); clf;
        set(gcf,'Color','White');
        clf(); fontsize = 20;
        h=contourf(znew,tnew,area_new,'edgecolor','none');
%         contourf(znew,tnew,area_new,'edgecolor','none');
        set(gca,'fontsize',fontsize','fontname','times');
%         set(gca,'XTick',linspace(newzmin,newzmax,6),'YTick',linspace(newtmin,newtmax,6));
%         set(gca,'YTickLabels',[0 50 650 700 750])

        % Load the cool-warm colormap
        cmap = load('CoolWarmFloat257.csv');
        colormap(cmap);
        cmin = min(area_new(:));
        cmax = max(area_new(:));
        caxis([cmin,cmax]);
        ch=colorbar();
        set(ch,'Ticks',linspace(cmin,cmax,5));




