% data_dir = 'H:\Numerics\KP\_tmax_125_Lx_300_Nx_512_Ly_200_Ny_256_bndry_condns_periodic_init_condns__solisegment__au_1_qu_0_ad_0_qd_0_x0_100_y0_0_w_100\';
% data_dir = 'H:\Numerics\KP\_tmax_150_Lx_800_Nx_512_Ly_400_Ny_256_bndry_condns_periodic_init_condns__solitest_true_soln\';
data_dir = 'H:\Numerics\KP\_tmax_150_Lx_800_Nx_1024_Ly_400_Ny_512_bndry_condns_periodic_init_condns__solitest_true_soln_tstart_20\';
num_on = 1;
save_on = 1;
%% ICs and Whitham theory predictions of evolution
    sam = 1; qm = 0;   x0 = 400; w = 100;
    sau = 0  ; qu =  sam;
    sad = 0  ; qd = -sam;
    up_lead_V = 2*qu - 2/3*sau;
    up_back_V = 2*qm - 2/3*sam;
    lo_back_V = 2*qm + 2/3*sam;
    lo_lead_V = 2*qd + 2/3*sad;
    % Initial condition of a RW started at time tstart
    tstart = 20;
%     Lx = 300; Nx = 2^8;
%     Ly = 200; Ny = Nx;
%     %% x and y vectors
%     x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
%     y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
%     [X,Y] = meshgrid(x,y);
load([data_dir,'parameters.mat'],'Nx','Ny','Lx','Ly','t');

%% x and y vectors
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [X,Y] = meshgrid(x,y);

if num_on
    
%% Process numerics
    %% Preallocate arrays
    up_lead = zeros(size(t));
    up_back = up_lead;
    lo_back = up_lead;
    lo_lead = up_lead;
    mean_u  = up_lead;
    max_u   = up_lead;
    ffty    = fftshift(y);
    %% denote amplitude cutoffs
    cutoff  = 0.05;
    kl = 0;
    kb = 0;
    li = length(t);
    bi = length(t);
	%% Extract leading/trailing edges, mean
    for ii = 1:length(t);
        load([data_dir,num2str(ii-1,'%05d.mat')],'u','tnow');
        u = gather(u);
        umax =  max(u,[],2);
        mean_u(ii) = mean(u(:));
        max_u(ii)  = max(u(:));
        try
            lo_back(ii) = y(find(abs(1-umax)<cutoff,1,'first'));
            up_back(ii) = y(find(abs(1-umax)<cutoff,1,'last'));
        catch
            if kb==0
                bi = ii;
            end
            kb = kb+1;
            continue;
        end
        try
            lo_lead(ii) = ffty(find(abs(  fftshift(umax))<cutoff,1,'last'));
            up_lead(ii) = ffty(find(abs(  fftshift(umax))<cutoff,1,'first'));
        catch
            if kl==0
                li = ii;
            end
            kl = kl+1;
            continue;
        end
    end
    %% Max useful time
        ti = 55;%round(length(t)/2);% min([bi,li,length(t)]);
   %% Plot characteristic speed comparison
    figure(1); clf;
    subplot(2,2,1);
        plot(t(1:ti),up_lead(1:ti),'.',...
             t(1:ti),up_lead_V*t(1:ti)+up_lead(1),'-');
        title('Upper leading edge');
	subplot(2,2,2);
        plot(t(1:ti),up_back(1:ti),'.',...
             t(1:ti),up_back_V*t(1:ti)+up_back(1),'-');
        title('Upper trailing edge');
	subplot(2,2,3);
        plot(t(1:ti),lo_back(1:ti),'.',...
             t(1:ti),lo_back_V*t(1:ti)+lo_back(1),'-');
        title('Lower trailing edge');
    subplot(2,2,4);
        plot(t(1:ti),lo_lead(1:ti),'.',...
             t(1:ti),lo_lead_V*t(1:ti)+lo_lead(1),'-');
        title('Lower leading edge');
	% Approximation to soliton maximum
    tasy =  t(t>75-tstart)/w;
    aasy = 1/4*(3./tasy).^(2/3) +...
            3^(1/3)/8*(1./tasy).^(4/3) +...
            3./(64*tasy.^2);
	figure(4); clf;
        subplot(2,1,1);
            plot(t,mean_u,'.');
        subplot(2,1,2);
            plot(t,max_u,'.',...
                tasy*w, aasy,'-');
end

return;
    %% Plot final useful time; compare to Whitham
    if num_on
%         x0 = soli.x0;
        [f2,ax2]=plot_nice_contour(data_dir,ti-1,2);
    else
        t = 50;
        [f2,ax2]=plot_nice_contour(data_dir,t,2);
%         [f2,ax2] = plot_nice_contour([pwd,filesep],t,2);
            ti = 1;
            title(['Time: ',num2str(t)]);
    end
    set(ax2,'XLim',[Lx/4 3*Lx/4],'YLim',[-Ly Ly],'CLim',[0 0.85]);

        set(f2,'Color','w'); set(gca,'fontsize',20);
        drawnow;
        if save_on
            print('seg_soli_numerics.png','-dpng');
        end

        
        [ uasy ] = vertical_segment(sam^2,sad^2,sau^2,...
                                     qm,qu,qd,...
                                     gather(x0),0,w,Lx,Ly,Nx,Ny,tstart+t(ti));
    U1 = uasy.us(X,Y);
    figure(3); clf; set(gcf,'Color','w');
    contourf(X,Y,U1,50,'edgecolor','none'); 
    set(gca,'CLim',get(ax2,'CLim'),'fontsize',20);
        set(gca,'XLim',get(ax2,'XLim'),'YLim',get(ax2,'YLim'));
        set(gcf,'Color','w'); set(gca,'fontsize',20);
        xlabel('x'); ylabel('y');
%     title(['Whitham theory, t=',num2str(t(ti))]);
    colormap(load('CoolWarmFloat257.csv')); colorbar;
    if save_on
            print('seg_soli_whitham.png','-dpng');
    end
