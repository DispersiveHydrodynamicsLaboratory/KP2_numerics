data_dir = '/Users/appm_admin/Documents/Numerics/KP/_tmax_150_Lx_800_Nx_1024_Ly_400_Ny_512_bndry_condns_periodic_init_condns__soliseg_check_tw_20/';
num_on = 1;
save_on = 0;
%% ICs and Whitham theory predictions of evolution
    sam = 1; qm = 0;   x0 = 200; w = 150;
    sau = 0  ; qu =  sam;
    sad = 0  ; qd = -sam;
    up_lead_V = 2*qu - 2/3*sau;
    up_back_V = 2*qm - 2/3*sam;
    lo_back_V = 2*qm + 2/3*sam;
    lo_lead_V = 2*qd + 2/3*sad;
    
    Lx = 800; Nx = 2^10;
    Ly = 400; Ny = 2^9;
    %% x and y vectors
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [X,Y] = meshgrid(x,y);

if num_on
    t  = 1:150;
%% x and y vectors
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [X,Y] = meshgrid(x,y);
    
%% Process numerics
    %% Preallocate arrays
    up_lead = zeros(size(t));
    up_back = up_lead;
    lo_back = up_lead;
    lo_lead = up_lead;
    mean_u  = up_lead;
    u_max   = up_lead;
    ffty    = fftshift(y);
    %% denote amplitude cutoffs
    cutoff  = 0.05;
    kl = 0;
    kb = 0;
    li = length(t);
    bi = length(t);
	%% Extract leading/trailing edges, mean
    for ii = 2:length(t);
        load([data_dir,num2str(ii-1,'%05d.mat')],'u');
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
%             continue;
        end
        try
            lo_lead(ii) = ffty(find(abs(  fftshift(umax))<cutoff,1,'last'));
            up_lead(ii) = ffty(find(abs(  fftshift(umax))<cutoff,1,'first'));
        catch
            if kl==0
                li = ii;
            end
            kl = kl+1;
%             continue;
        end
        u_max(ii) = findpeaks(u(Ny/2,:),x,'SortStr','descend','Npeaks',1);
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
% 	figure(4); clf;
%         plot(t,mean_u,'.');
end

t = 50;
   load([data_dir,num2str(t,'%05d.mat')],'u');
   figure(1); clf;
    plot(x,u(Ny/2,:));



    %% Plot final useful time; compare to Whitham
    if num_on
%         x0 = soli.x0;
        [f2,ax2]=plot_nice_contour(data_dir,ti-1,2);
    else
        t = 10;
        [f2,ax2] = plot_nice_contour([data_dir,filesep],t,2);
            ti = 1;
        
    end
%     set(ax2,'XLim',[0 Lx/2],'YLim',[-Ly Ly],'CLim',[0 0.85]);
        set(f2,'Color','w'); set(gca,'fontsize',20);
        title('');
        if save_on
            print('seg_soli_numerics.png','-dpng');
        end

    figure(3); clf; set(gcf,'Color','w');
    a1 = @(x,y,t) 9/64 * ( 2- y/t).^2 .*( (-2/3*t) <= y ).*( y < ( 2*t) ) + ...
                 1                    .*( (-2/3*t) > y );
    a2 = @(x,y,t) 9/64 * ( 2+ y/t).^2 .*( ( 2/3*t) >= y ).*( y > (-2*t) )+ ...
                 1                    .*( ( 2/3*t) < y );
    q1 = @(x,y,t) ( 1-sqrt(a1(x,y,t)));
	q2 = @(x,y,t) (-1+sqrt(a2(x,y,t)));
    
    a = @(x,y,t) a1(x,y-w/2,t).*(y>=0) + a2(x,y+w/2,t).*(y<0);
    q = @(x,y,t) q1(x,y-w/2,t).*(y>=0) + q2(x,y+w/2,t).*(y<0);
    qy = q1((X-x0),(Y-w/2),t(ti)).*(Y-w/2) + q2((X-x0),(Y+w/2),t(ti)).*(Y<0).*(Y+w/2);
    
    U1 = a((X-x0),Y,t(ti)).*sech( sqrt(a((X-x0),Y,t(ti))/12).*...
          ((X-x0) + qy - (a((X-x0),Y,t(ti))/3 + q((X-x0),Y,t(ti)).^2)*t(ti))).^2;
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

    
    