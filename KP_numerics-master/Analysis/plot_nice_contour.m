function [ fh,ax ] = plot_nice_contour(data_dir,tind, fignum )
%PLOT_NICE_CONTOUR Summary of this function goes here
%   Detailed explanation goes here
cnum=50;

load([data_dir,'parameters.mat'],'Nx','Ny','Lx','Ly');
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
load([data_dir,num2str(tind,'%05d'),'.mat'],'u','tnow');
fh = figure(fignum); clf;
    contourf(x,y,u,cnum,'edgecolor','none')
    ax = gca;
    colormap(load('CoolWarmFloat257.csv'));
    colorbar;
    xlabel('x'); ylabel('y'); title(['t: ',num2str(tnow)]);
end

