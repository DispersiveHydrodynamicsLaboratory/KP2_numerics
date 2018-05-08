% data_dir = 'H:\Numerics\KP\_tmax_125_Lx_300_Nx_512_Ly_200_Ny_256_bndry_condns_periodic_init_condns__solisegment__au_1_qu_0_ad_0_qd_0_x0_100_y0_0_w_100\';
% data_dir = 'H:\Numerics\KP\_tmax_150_Lx_800_Nx_512_Ly_400_Ny_256_bndry_condns_periodic_init_condns__solitest_true_soln\';
% data_dir = 'H:\Numerics\KP\_tmax_150_Lx_800_Nx_1024_Ly_400_Ny_512_bndry_condns_periodic_init_condns__solitest_true_soln_tstart_20\';
% rip;
data_dir = '/Volumes/Data Storage/Numerics/KP/_tmax_150_Lx_800_Nx_1024_Ly_400_Ny_512_bndry_condns_periodic_init_condns__solitest_true_soln/';
num_on = 1;
save_on = 0;
fs = 16; %fontsize for nice pictures
%% ICs and Whitham theory predictions of evolution
    sam = 1; qm = 0;   x0 = 200; w = 100; thodo = 3/4*w;
    sau = 0  ; qu =  sam;
    sad = 0  ; qd = -sam;
    up_lead_V = 2*qu - 2/3*sau;
    up_back_V = 2*qm - 2/3*sam;
    lo_back_V = 2*qm + 2/3*sam;
    lo_lead_V = 2*qd + 2/3*sad;
    % Initial condition of a RW started at time tstart
    tstart = 30; 
    % True time for transition to hodograph transf
    tth = thodo-tstart;
%     Lx = 300; Nx = 2^8;
%     Ly = 200; Ny = Nx;
%     %% x and y vectors
%     x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
%     y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
%     [X,Y] = meshgrid(x,y);
load([data_dir,'parameters.mat'],'Nx','Ny','Lx','Ly','t','soli');
%% x and y vectors
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [X,Y] = meshgrid(x,y);
%% Generate text file with parameters
fileID = fopen('seg_soli_params.txt','w');
    fprintf(fileID,['Dir Name: ',data_dir,'\n']);
    fprintf(fileID,['Soliton Amplitudes:   ',...
        'Top: ',num2str(sau^2),'   Mid: ',num2str(sam^2),'   Bot: ',num2str(sad^2),'\n']);
    fprintf(fileID,['Profile start: ',num2str(tstart), '   Middle portion width: ',num2str(w)]);
    fclose(fileID);
    
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
    cutoff_lead  = 0.05;
    cutoff_back  = 0.1;
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
            lo_back(ii) = y(find(abs(1-umax)<cutoff_back,1,'first'));
            up_back(ii) = y(find(abs(1-umax)<cutoff_back,1,'last'));
        catch
            if kb==0
                bi = ii;
            end
            kb = kb+1;
%             continue;
        end
        try
            lo_lead(ii) = ffty(find(abs(  fftshift(umax))<cutoff_lead,1,'last'));
            up_lead(ii) = ffty(find(abs(  fftshift(umax))<cutoff_lead,1,'first'));
        catch
            if kl==0
                li = ii;
            end
            kl = kl+1;
%             continue;
        end
    end
    %% Max useful time
    [~,ti_th] = min(abs(t-(thodo-tstart)));
        ti = min([bi,li,ti_th]);
    %% Plot characteristic speed comparison
    figure(1); clf;
    set(gcf,'Color','White');
    sax(1)=subplot(2,2,1);
        plot(t(1:ti),up_lead(1:ti),'.',...
             t(1:ti),up_lead_V*t(1:ti)+up_lead(1),'-');
        title('Upper leading edge');
        xlabel('t'); ylabel('y');
	sax(2)=subplot(2,2,2);
        plot(t(1:ti),up_back(1:ti),'.',...
             t(1:ti),up_back_V*t(1:ti)+up_back(1),'-');
        title('Upper trailing edge');
        xlabel('t'); ylabel('y');
	sax(3)=subplot(2,2,3);
        plot(t(1:ti),lo_back(1:ti),'.',...
             t(1:ti),lo_back_V*t(1:ti)+lo_back(1),'-');
        title('Lower trailing edge');
        xlabel('t'); ylabel('y');
    sax(4)=subplot(2,2,4);
        plot(t(1:ti),lo_lead(1:ti),'.',...
             t(1:ti),lo_lead_V*t(1:ti)+lo_lead(1),'-');
        title('Lower leading edge');
        xlabel('t'); ylabel('y');
    set(sax,'Fontsize',fs);
    if save_on
        print('seg_soli_speeds.png','-dpng');
    end
    %% Approximation to soliton maximum
    tasy =  (t(t>=tth)+tstart)/w;
    aasy = 1/4*(3./tasy).^(2/3) +...
            3^(1/3)/8*(1./tasy).^(4/3) +...
            3./(64*tasy.^2);
	figure(2); clf;
            plot(t,max_u,'.',...
                (tasy*w)-tstart, aasy,'-');
            hold on; plot([tth tth],get(gca,'YLim'));
            xlabel('t'); ylabel('a(t,y=0)');
            title(['Soliton Amplitude Decay']);
            legend('Numerics','Hodograph Approximation','Transition: Simple Wave to Hodograph');
            set(gcf,'Color','White');
            set(gca,'Fontsize',fs);
    if save_on
        print('seg_soli_ampl.png','-dpng');
    end
    
    %% Check for zero mean, conserved quantities (WIP)
	figure(3); clf;
        plot(t,mean_u,'.');
        xlabel('t'); ylabel('mean background');
        title(['Maximum mean: ',num2str(max(mean_u(:)))]);
        set(gcf,'Color','White');
            set(gca,'Fontsize',fs);
    if save_on
        print('seg_soli_conserved_quants.png','-dpng');
    end
end

    %% Plot final useful time; compare to Whitham
    if num_on
        [f4,ax4]=plot_nice_contour(data_dir,ti-1,4);
        title(['Time: ',num2str(t(ti)),' Numerics']);
    else
        t = 50;
        [f4,ax4]=plot_nice_contour(data_dir,t,4);
%         [f2,ax2] = plot_nice_contour([pwd,filesep],t,2);
            ti = 1;
            title(['Time: ',num2str(t),' Numerics']);
    end
    set(ax4,'XLim',[x0-100 x0+100],'YLim',[-Ly Ly],'CLim',[0 1]);
        set(f4,'Color','w'); set(gca,'fontsize',fs);
        drawnow;
        if save_on
            print('seg_soli_numerics.png','-dpng');
        end

        
        [ uasy ] = vertical_segment(sam^2,sad^2,sau^2,...
                                     qm,qu,qd,...
                                     gather(x0),0,w,Lx,Ly,Nx,Ny,tstart+t(ti));
    U1 = uasy.us(X,Y);
    figure(5); clf; set(gcf,'Color','w');
    contourf(X,Y,U1,50,'edgecolor','none'); 
    set(gca,'CLim',get(ax4,'CLim'),'fontsize',fs);
        set(gca,'XLim',get(ax4,'XLim'),'YLim',get(ax4,'YLim'));
        set(gcf,'Color','w'); set(gca,'fontsize',fs);
        xlabel('x'); ylabel('y');
        title(['Time: ',num2str(t(ti)),' Whitham']);
    colormap(load('CoolWarmFloat257.csv')); colorbar;
    if save_on
            print('seg_soli_whitham.png','-dpng');
    end
