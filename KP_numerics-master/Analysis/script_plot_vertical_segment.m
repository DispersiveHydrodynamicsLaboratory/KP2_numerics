%% Master figure generator for KP2 line soliton segment
close all;
% Data directory 
if strcmp(getenv('computername'),'CAKELIE2012')
    main_dir = 'H:\Numerics\KP\';
else
    main_dir = pwd;
end
run_dir  = '_tmax_150_Lx_800_Nx_1024_Ly_400_Ny_512_bndry_condns_periodic_init_condns__solitest_true_soln\';
data_dir = [main_dir, run_dir];

% Settings
save_on = 1; % set to 1 to save plots generated, 0 to now
lw = 1;      % line width
ms = 3;      % markersize
fstem = 'seg_soli_'; % Stem for naming convention

%% ICs and Whitham theory predictions of evolution (hard-coded for ease of use)
    sam = 1  ; qm = 0;  
    sau = 0  ; qu =  sam;
    sad = 0  ; qd = -sam;
    up_lead_V = 2*qu - 2/3*sau;
    up_back_V = 2*qm - 2/3*sam;
    lo_back_V = 2*qm + 2/3*sam;
    lo_lead_V = 2*qd + 2/3*sad;
    % Center of initial soliton, soliton width, Whitham prediction for
    % transition
    x0 = 200; y0=0; w = 100; thodo = 3/4*w;
    % Initial condition of a RW started at time tstart
    tstart = 30; 
    % True time for transition to hodograph transf
    tth = thodo-tstart;
    
%% x and y vectors
load([data_dir,'parameters.mat'],'Nx','Ny','Lx','Ly','t');
    x = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    y = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [X,Y] = meshgrid(x,y);
%% Representative time for simple wave theory (hardcoded)
    t_useful = 40;
    [~,ti] = min(abs(t-t_useful));

%% Generate text file with numerics information
fileID = fopen('seg_soli_params.txt','w');
    fprintf(fileID,['Dir Name: ',run_dir,'\n']);
    fprintf(fileID,['Soliton Amplitudes:   ',...
        'Top: ',num2str(sau^2),'   Mid: ',num2str(sam^2),'   Bot: ',num2str(sad^2),'\n']);
    fprintf(fileID,['Soliton Angles of Propagation:   ',...
        'Top: ',num2str(qu),'   Mid: ',num2str(qm),'   Bot: ',num2str(qd),'\n']);
    fprintf(fileID,['Profile start: ',num2str(tstart), '   Middle portion height: ',num2str(w),'\n']);
    fprintf(fileID,['Time of transition to hodograph: ',num2str(tth),'\n']);
    fclose(fileID);

%% Process numerics
    % Preallocate arrays
    up_lead = zeros(size(t));
    up_back = up_lead;
    lo_back = up_lead;
    lo_lead = up_lead;
    mean_u  = up_lead;
    max_u   = up_lead;
    dx_u    = up_lead;
    ffty    = fftshift(y);
    % denote amplitude cutoffs
    cutoff_lead  = 0.05;
    cutoff_back  = 0.1;
    fileID = fopen('seg_soli_params.txt','a');
        fprintf(fileID,['Leading Edge Determination: find where soliton ampliude is: ',num2str(cutoff_lead),...
                    ', where leading edges are edges moving away from the center. \n']);
    	fprintf(fileID,['Trailing Edge Determination: find where soliton ampliude is: ',num2str(cutoff_back),...
                    ', where trailing edges are edges moving towards the center. \n']);
	fclose(fileID);
    kl = 0;
    kb = 0;
    li = length(t);
    bi = length(t);
	% Extract leading/trailing edges, mean, max, area integral of u_x
    for ii = 1:length(t);
        load([data_dir,num2str(ii-1,'%05d.mat')],'u','tnow');
        u = gather(u);
        umax =  max(u,[],2);
        mean_u(ii) = mean(u(:));
        max_u(ii)  = max(u(:));
        % u_x calculation via second-order centered FD (periodic BCs)
        ux = ( u(:,[2:end,1]) - u(:,[end,1:end-1]) )/(2*(2*Lx/Nx));
        dx_u(ii) = trapz(y,trapz(x,ux,2));
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
%% Fig. 1: Plot characteristic speed comparison
inds = 1:3:(ti-1); % Too many points looks cluttered
fspeed = figure(1); clf;
set(gcf,'Color','White');
pspeed = plot(t(inds),up_lead_V*t(inds)+up_lead(1),'k-',...% Leading ref line
              t(inds),up_back_V*t(inds)+up_back(1),'k-',...% Trailing ref line
              t(inds),up_lead(inds),'b^',...     %Upper leading edge;
              t(inds),up_back(inds),'b^',...     %Upper trailing edge;
              t(inds),abs(lo_back(inds)),'rs',...%Lower trailing edge;
              t(inds),abs(lo_lead(inds)),'rs',...%Lower leading edge;
              'LineWidth',lw,'MarkerSize',ms); 
     for ii = 3:4
         pspeed(ii).MarkerFaceColor = 'b';
         pspeed(ii+2).MarkerFaceColor = 'r';
     end
    legend(pspeed([1 3 5]),...
           {'Whitham','Upper Half Plane','Lower Half Plane'},...
           'Location','NorthWest');
    leg_on = 1; ch_on = 0;
    xlabel('t'); ylabel('y');
fspeed = make_num_plots_for_pub(fspeed,[fstem,'speeds'],ch_on,leg_on,save_on);

%% Fig. 2: Approximation to soliton maximum; Comparison to Hodograph
hodo_a  = @(tasy) 1/4*(3./tasy).^(2/3) +...
        3^(1/3)/8*(1./tasy).^(4/3) +...
        3./(64*tasy.^2);
tasy_pure  =  (t(t>=tth)+tstart)/w;
% Curve Fitting gets an estimated phase shift of \approx 15
tshift = 15;
tasy_shift =  (t(t>=tth-tshift)+tstart+tshift)/w;
fhodot = figure(2); clf;
        plot(t,max_u,'k.',...
            (tasy_pure*w)-tstart, hodo_a(tasy_pure),'b-',...
            (tasy_shift*w)-(tstart+tshift), hodo_a(tasy_shift),'r-',...
            'LineWidth',lw,'MarkerSize',ms);
        xlabel('t'); ylabel('a(t,y=0)');
        legend('Numerics','Hodograph','Shifted Hodograph','Location','NorthEast');
        leg_on = 1; ch_on = 0;
fhodot = make_num_plots_for_pub(fhodot,[fstem,'hodo'],ch_on,leg_on,save_on);
    
%% Fig. 3: Check for zero mean (print to file), u_x conserved throughout run (figure)
fcons = figure(3); clf;
    fileID = fopen('seg_soli_params.txt','a');
        fprintf(fileID,['Maximum mean: ',num2str(max(mean_u(:))),'\n']);
        fclose(fileID);
        plot(t,dx_u,'.');
        xlabel('t'); ylabel('$u_x$');
    set(gcf,'Color','White');
        leg_on = 0; ch_on = 0;
fcons = make_num_plots_for_pub(fcons,[fstem,'cons'],ch_on,leg_on,save_on);

%% Fig. 4: Plot representative time: Numerics
ti = 41; % For most trials, t(ti) = ti-1
fileID = fopen('seg_soli_params.txt','a');
    fprintf(fileID,['Contour plot time: ',num2str(t(ti))]);
fclose(fileID);
cnum=50;
load([data_dir,num2str(ti,'%05d'),'.mat'],'u','tnow');
fnum = figure(4); clf;
    contourf(x,y,u,cnum,'edgecolor','none');
	% Compare to Whitham Theory
    [ uasy ] = vertical_segment(sam^2,sad^2,sau^2,...
                                 qm,qu,qd,...
                                 x0,y0,w,Lx,Ly,Nx,Ny,tstart+t(ti));
    % Overlay the Whitham theory prediction with theta=0
    [~,theta0] = min(abs(uasy.th.f(X-x0,Y-y0)),[],2);
    hold on;
    plot(x(theta0),y,'k--','LineWidth',lw);

    axnum = gca;
    colormap(load('CoolWarmFloat257.csv'));
    colorbar;
    xlabel('x'); ylabel('y');
    set(axnum,'XLim',[x0-100 x0+100],'YLim',[-Ly Ly],'CLim',[0 1]);
        leg_on = 0; ch_on = 1;
fnum = make_num_plots_for_pub(fnum,[fstem,'numerics'],ch_on,leg_on,save_on);

%% Fig. 5: Plot representative time: Whitham
[ uasy ] = vertical_segment(sam^2,sad^2,sau^2,...
                             qm,qu,qd,...
                             x0,y0,w,Lx,Ly,Nx,Ny,tstart+t(ti));
U1 = uasy.us(X,Y);
fwhit = figure(5); clf; set(gcf,'Color','w');
    contourf(X,Y,U1,50,'edgecolor','none'); 
        set(gca,'XLim',get(axnum,'XLim'),'YLim',get(axnum,'YLim'));
        set(gca,'CLim',get(axnum,'CLim'));
        set(gcf,'Color','w'); 
        xlabel('x'); ylabel('y');
    colormap(load('CoolWarmFloat257.csv')); colorbar; drawnow;
        leg_on = 0; ch_on = 1;
fwhit = make_num_plots_for_pub(fwhit,[fstem,'whitham'],ch_on,leg_on,save_on);

%% Fig. 6: Compare soliton amplitudes for representative time
a_num  = max(u,[],2);
a_whit = max(U1,[],2);
    fampl = figure(6); clf;
    set(gcf,'Color','White');
    plot(y,a_num,...
         y,a_whit,...
         'LineWidth',lw,'MarkerSize',ms);
		legend({'Numerics','Whitham'},'Location','NorthEast');
        leg_on = 1; ch_on = 0;
        xlabel('t'); ylabel('y');
    fampl = make_num_plots_for_pub(fampl,[fstem,'ampls'],ch_on,leg_on,save_on);



