function plot_data_fun_2D(loaddir,varargin)
load([loaddir,'parameters.mat'],'t','Nx','Ny','dx','dy','xmax','ymax'); %,'f');
  Lx = xmax; Ly = ymax;
    x = (2*Lx/Nx)*(-Nx/2:Nx/2-1)';
    y = (2*Ly/Ny)*(-Ny/2:Ny/2-1)';
if nargin>1
    nplts = varargin{1};
    if nargin>2
        tmax = varargin{2};
        if nargin > 3
            tmin = varargin{3};
        else
            tmin = 0;
        end
    else
        tmax = Inf;
        tmin = 0;
    end
else
    nplts = 6;
    tmax  = Inf;
    tmin  = 0;
end

% Find maximum t index
tind = length(t)-1;
for ii=1:length(t)-1
  [fid,foo] = fopen(strcat(loaddir,num2str(ii,'%05d'),'.mat'),'r');
  if fid == -1 % File does not exist
    tind = ii-1;
    disp(['Maximum time = ',num2str(t(tind))]);
    break;
  end
  fclose(fid);
end

cmap = load('CoolWarmFloat257.csv');
tmaxind = find(t<=tmax,1,'last');
tminind = find(t>=tmin,1,'first');

t = t(1:min(tind,tmaxind));

% Plotting stuff
plot_figs = [1   % Time sequence of n plots
             0   % Boundary deviations
             1   % Initial vs final time
             0]; % Conserved quantity


fontsize = 12;
% Extract domain over which to plot
[foo,xmini] = min(abs(x + xmax));
[foo,xmaxi] = min(abs(x - xmax));
xp = xmini:xmaxi;
[foo,ymini] = min(abs(y + ymax));
[foo,ymaxi] = min(abs(y - ymax));
yp = ymini:ymaxi;


if plot_figs(1)
    %toutind = 1:6;
    if nplts ~= 1
        toutind = round(linspace(max(1,tminind),length(t)-1,nplts-1)); %linspace(round((length(t)-1)/(nplts-1))
    else
        toutind = tminind;
    end
    % Find max and min of solution to plot
%     if exist('f','var')
%         umax = max(f(x(xp)));
%         umin = min(f(x(xp)));
%     else
        load(strcat(loaddir,num2str(0,'%05d'),'.mat'),'u_init');
        umax = max(u_init(:)); %max(f(z(zp)));
        umin = min(u_init(:)); %min(f(z(zp)));
%     end
    for ii=1:max(nplts-1,1)
        load(strcat(loaddir,num2str(toutind(ii),'%05d'),'.mat'),'u','tnow','inc');
        if max(u(:)) > umax
            umax = max(u(:));
        end
        if min(u(:)) < umin
            umin = min(u(:));
        end
    end
    
    figure(4)
    clf()
    % Plot initial condition
    load(strcat(loaddir,num2str(0,'%05d'),'.mat'),'u_init');
    h(1)=subplot(nplts,1,1);
        contourf(x(xp),y(yp),u_init(yp,xp),100,'edgecolor','none');
    set(gca,'fontsize',fontsize,'fontname','times');
    ylabel('$u$','interpreter','latex');
    axis([x(xmini),x(xmaxi),y(ymini),y(ymaxi)]);
    colormap(cmap);
    caxis([umin umax])
    title({['$t = 0$']},'interpreter','latex');
    for ii=1:max(nplts-1,1)
        load(strcat(loaddir,num2str(toutind(ii),'%05d'),'.mat'),'u','tnow','inc');
        if nplts ~= 1
            h(ii+1) = subplot(nplts,1,ii+1);
        else
            h(1) = subplot(1,1,1);
        end
        contourf(x(xp),y(yp),real(u(yp,xp)),100,'edgecolor','none');
        caxis([umin umax])
        set(gca,'fontsize',fontsize,'fontname','times');
        ylabel('$y$','interpreter','latex');
        if ii == nplts-1
            xlabel('$x$','interpreter','latex');
%             colorbar; 
        end
        axis([x(xmini),x(xmaxi),y(ymini),y(ymaxi)]);
        if nplts ~= 1
            title({['$t = ',num2str(t(toutind(ii)+1)),'$']},'interpreter','latex');
        else
            title({['$t = ',num2str(t(1)),'$']},'interpreter','latex');
        end
    end
end
linkaxes(h,'xy')

if plot_figs(2)
    % Extract solution values adjacent to boundaries
    A1 = zeros(1,tind);
    AN = zeros(1,tind);
    tp = zeros(1,tind);
    for ii=1:tind
        load(strcat(loaddir,num2str(ii,'%05d'),'.mat'),'A','tnow','inc');
        A1(ii) = A(1);
        AN(ii) = A(Nz);
        tp(ii) = tnow;
    end
    
    figure(4)
    clf()
    subplot(1,2,1);
    plot(tp,g0(tp),'b-',...
         tp,A1,'r--');
    set(gca,'fontsize',fontsize,'fontname','times');
    legend({['$A(0,t)$'],['$A(z_1,t)$']},'interpreter','latex',...
           'location','northwest');
    xlabel('$t$','interpreter','latex');
    subplot(1,2,2);
    plot(tp,g1(tp),'b-',...
         tp,AN,'r--');
    set(gca,'fontsize',fontsize,'fontname','times');
    legend({['$A(L,t)$'],['$A(z_N,t)$']},'interpreter','latex',...
           'location','northwest');
    xlabel('$t$','interpreter','latex');
end

% Initial and final times
if plot_figs(2)
    load(strcat(loaddir,num2str(tind,'%05d'),'.mat'),'A','tnow','inc');
    
    figure(4)
    clf()
    subplot(1,2,1);
    plot(x,f(x),'b-',...
         x,A,'r--');
    set(gca,'fontsize',fontsize,'fontname','times');
    legend({['$A(z,0)$'],['$A(z,t_{\rm f})$']},'interpreter','latex',...
           'location','northwest');
    xlabel('$z$','interpreter','latex');
    subplot(1,2,2);
    plot(x,f(x)'-A,'b-');
    set(gca,'fontsize',fontsize,'fontname','times');
    xlabel('$z$','interpreter','latex');
    ylabel('$A(z,0)-A(z,t_{\rm f})$','interpreter','latex');
end


