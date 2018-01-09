function[fh] = make_contours_for_pub(u,z,t,tmin,tmax,...
                            zmin,zmax,cmin,cred,cmax,...
                            fignum,fhname,slant)
% TODO: many of these could be inputs
debug_on = 0;
save_on  = 1;
Ninterp  = 500;
outerpos = 0; 
fontsize = 16;
fontname = 'helvetica';
w1       = 8*2.54;
h1       = 8*2.54;%w1/1.618;
units    = 'centimeters';
fspec    = '-dpdf';
res      = '-r0';
xlabelsize = 6;
ylabelsize = 6;
% Rename space and time vectors for consistency
    T    = t;
    zold = z;

% Set up grid to interpolate TO
    zmin = max(min(zold),zmin);
    zmax = min(max(zold),zmax);
    tmin = max(min(T),tmin);
    tmax = min(max(T),tmax);
    znew = linspace(zmin, zmax, Ninterp);
    tnew = linspace(tmin, tmax, Ninterp);
    
%   Interpolate to new vectors
    unew = interp2(zold,T,u,znew,tnew','spline');

disp('plotting...');
        % DALTON
        disp('Plotting contour plot...');
        fh = figure(fignum);clf
        if cred
           cv = [linspace(min(unew(:)),cred,25) linspace(cred*1.01,max(unew(:)),25)];
        else
           cv = 50;
        end
            contourf(znew,tnew,unew,cv,'edgecolor','none'); 
            % Add line showing trajectory
            lwidth = 0.5;
            hold on;
                ah = gca;
                set(ah,'fontsize',fontsize,'fontname',fontname);
                set(ah,'XTick',floor(linspace(zmin,zmax,xlabelsize)/10)*10);
                set(ah,'TickLabelinterpreter','latex');
                set(ah,'YTick',floor(linspace(tmin,tmax,ylabelsize)/10)*10);
            cmin = max(min(unew(:)),cmin);
            cmax = min(max(unew(:)),cmax);
            cmap = load('CoolWarmFloat257.csv');
            % Edit colormap for better visibility
            if cred
                cmap = cmap_edit(cred,cmin,cmax,cmap);
            end
            
            colormap(cmap);
            caxis([cmin,cmax]);
            ch=colorbar('Fontsize',fontsize,'Fontname',fontname);
            set(ch,'TickLabelInterpreter','latex');
            
            xlabel('x','HorizontalAlignment','center',...
                            'Interpreter','latex','Fontsize',fontsize);
            yname = ylabel('y',...
                            'Interpreter','latex','Fontsize',fontsize);

        % ---- change figure window ---- 
        set(fh,'units',units,'color','white')
        pos_fh = get(fh,'position'); % figure window
        set(fh,'position',[1 0,w1,h1])
        drawnow;
        % ---- change plot axes ----
        axs = 0.6; %normalized to width
        set(ah,'units',units,'fontsize',fontsize)
        pos_ah = get(ah,'position'); % plot axes
        set(ah,'position',[pos_ah(1:2),axs*w1,pos_ah(4)])
        % ---- change colorbar axes -----
        set(ch,'units',units,'fontsize',fontsize)
        pos_ch = get(ch,'position'); % colorbar
        pos_ch(1) = pos_ah(1) + (axs*1.1)*w1;
        set(ch,'position',[pos_ch(1:2),0.5*pos_ch(3),pos_ch(4)])
        % ---- change ylabel position
        ynamepos = get(yname,'Position');
        set(yname,'Position',[(1-2/32)*ynamepos(1), ynamepos(2:3)]);
        % ---- adjust paper size for printing ----
        set(fh,'PaperUnits','centimeters',...
               'PaperPosition',[0,0,w1,h1],...
               'PaperSize',[w1,h1])
           
        % Print figures
        if save_on
            disp('Saving...');
                print(fh,fhname,fspec,'-r300');
                disp(['Contour plot saved to ',fhname,' as a ',fspec(3:end)]);
        end
        hold off;


