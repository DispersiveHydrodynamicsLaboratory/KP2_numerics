data_dir = 'H:\Numerics\KP\_tmax_125_Lx_300_Nx_512_Ly_200_Ny_256_bndry_condns_periodic_init_condns__solisegment__au_1_qu_0_ad_0_qd_0_x0_100_y0_0_w_100\';

load([data_dir,'parameters.mat']);
%% x and y vectors
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [X,Y] = meshgrid(x,y);
%% ICs and Whitham theory predictions of evolution
    sam = 1; qm = 0;
    sau = 0  ; qu = sam;
    sad = 0  ; qd = -sam;
    up_lead_V = 2*qu - 2/3*sau;
    up_back_V = 2*qm - 2/3*sam;
    lo_back_V = 2*qm + 2/3*sam;
    lo_lead_V = 2*qd + 2/3*sad;
    
%% Process numerics
    %% Preallocate arrays
    up_lead = zeros(size(t));
    up_back = up_lead;
    lo_back = up_lead;
    lo_lead = up_lead;
    mean_u  = up_lead;
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
        umax =  max(u,[],2);
        mean_u(ii) = mean(u(:));
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
        ti = min([bi,li,length(t)]);
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
	figure(4); clf;
        plot(t,mean_u,'.');
    %% Plot final useful time; compare to Whitham
    x0 = soli.x0;
    [f2,ax2]=plot_nice_contour(data_dir,ti-1,2);
    figure(3); clf;
    a = @(x,y,t) 9/64*(2-(y-50)/t).^2 .*(-2/3*t<=(y-50)).*((y-50)< 2*t) + ...
                 9/64*(2+(y+50)/t).^2 .*( 2/3*t>=y+50).*(y+50>-2*t) + ...
                 1                .*(-2/3*t>y-50).*(2/3*t<y+50);
             
    q = @(x,y,t) -1                 .*(-2*t>=y+50) +...
                 +1                 .*( 2*t<=y-50)+ ...
                 (1-3/8*(2-(y-50)/t)).*(-2/3*t<=y-50).*(y-50< 2*t) + ...
                 (-1+3/8*(2+(y+50)/t)).*( 2/3*t>=y+50).*(y+50>-2*t);
%     ti=75
    U = a((X-x0),Y,t(ti)).*sech( sqrt(a((X-x0),Y,t(ti))/12).*...
         ((X-x0) + q((X-x0),Y,t(ti)).*Y - (a((X-x0),Y,t(ti))/3 + q((X-x0),Y,t(ti)).^2)*t(ti))).^2;
    contourf(X,Y,U,50,'edgecolor','none'); 
    set(gca,'CLim',get(ax2,'CLim'));
    colormap(load('CoolWarmFloat257.csv')); colorbar;
     
