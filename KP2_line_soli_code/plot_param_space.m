%% Code for generating KP boundaries
save_on = 1;
fhname = 'KP2_param_x';
w1 = 8;
h1 = 8;
units = 'inches';
Lmax = 5;
fs = 36;
fontsize = fs;
x = -Lmax:0.01:Lmax;
y = x;
[X,Y] = meshgrid(x,y);

%% x-direction
fh = figure(1); clf; set(gcf,'Color','White');
    A = -(0.5*(X<Y) + (X>-Y)) ;
    contourf(X,Y,A);
    ah = gca;
    set(ah,'Fontname','helvetica','fontsize',fs);
    colormap(load('CoolWarmFloat257.csv'));
    t1 = text(Lmax/2,0,'1-RW, 2-shock','fontsize',fs/2,'HorizontalAlignment','center');
    t2 = text(-Lmax/2,0,'2-RW, 1-shock','fontsize',fs/2,'HorizontalAlignment','center');
    t3 = text(0,Lmax/2,'all-RW','Color','w','fontsize',fs/2,'HorizontalAlignment','center');
    t4 = text(0,-Lmax/2,'all-shock','Color','w','fontsize',fs/2,'HorizontalAlignment','center');
%     xlabel('a_l - a_r');
%     ylabel('q_l - q_r');
    axis square
    drawnow;
    
        % ---- change figure window ---- 
        set(fh,'units',units,'color','white')
        pos_fh = get(fh,'position'); % figure window
        set(fh,'position',[pos_fh(1:2),w1,h1])
        drawnow;
        % ---- change plot axes ----
        axs = 0.8; %normalized to width
        set(ah,'units',units,'fontsize',fontsize)
        pos_ah = get(ah,'position'); % plot axes
        set(ah,'position',[pos_ah(1:2),axs*w1,pos_ah(4)])
%         % ---- change colorbar axes -----
%         set(ch,'units',units,'fontsize',fontsize)
%         pos_ch = get(ch,'position'); % colorbar
%         pos_ch(1) = pos_ah(1) + (axs*1.1)*w1;
%         set(ch,'position',[pos_ch(1:2),0.5*pos_ch(3),pos_ch(4)])
%         % ---- change ylabel position
%         ynamepos = get(yname,'Position');
%         set(yname,'Position',[(1-2/32)*ynamepos(1), ynamepos(2:3)]);
        % ---- adjust paper size for printing ----
        set(fh,'PaperUnits',units,...
               'PaperPosition',[0,0,w1,h1],...
               'PaperSize',[w1,h1])
           
        % Print figures
        if save_on
            disp('Saving...');
                print(fh,fhname,'-dpdf');
                disp(['Contour plot saved to ',fhname]);
        end
        hold off;
