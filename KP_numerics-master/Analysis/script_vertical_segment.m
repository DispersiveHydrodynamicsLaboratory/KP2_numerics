%% ICs and Whitham theory predictions of evolution
    sam = 1; qm = 0;   x0 = 0; w = 100;
    sau = 0  ; qu =  sam;
    sad = 0  ; qd = -sam;
    % time
    tstart = 75;
    Nx = 2^8;
    Ny = 2^7;
    Lx = 400;
    Ly = 200;
%% x and y vectors
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [X,Y] = meshgrid(x,y);
    
        [ uasy ] = vertical_segment(1,0,0,...
                                     0,-1,1,...
                                     0,0,w,Lx,Ly,Nx,Ny,tstart);
    U1 = uasy.us(X,Y);
    figure(3); clf; set(gcf,'Color','w');
    contourf(X,Y,U1,50,'edgecolor','none'); 
%     set(gca,'CLim',get(ax2,'CLim'),'fontsize',20);
%         set(gca,'XLim',get(ax2,'XLim'),'YLim',get(ax2,'YLim'));
        set(gcf,'Color','w'); set(gca,'fontsize',20);
        xlabel('x'); ylabel('y');
%     title(['Whitham theory, t=',num2str(t(ti))]);
    colormap(load('CoolWarmFloat257.csv')); colorbar;
figure(2);contourf(X,Y,uasy.th.f(X-uasy.x0,Y- uasy.y0),50,'edgecolor','none'); caxis([-100 100]); colorbar;
hold on;
plot([x0 x0],[min(Y(:)) max(Y(:))],'k-');