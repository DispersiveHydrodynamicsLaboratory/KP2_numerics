function[fh] = make_num_plots_for_pub(fh,fhname,ch_on,leg_on,save_on)
% This function formats plots

% TODO: many of these could be inputs
fontsize = 7;
fontname = 'times';
w1       = 3*2.54;
h1       = 2*2.54;%w1/1.618;
units    = 'centimeters';
fspec    = '-dpdf';
res      = '-r0';

        figure(fh); set(fh,'Color','White');
        
        %% Set Axes 
            ah = gca;
            set(ah,'fontsize',fontsize,'fontname',fontname,...
                'TickLabelInterpreter','latex');
            if ch_on
                ch=colorbar('Fontsize',fontsize,'Fontname',fontname);
            end
            xlabel(ah.XLabel.String,'HorizontalAlignment','center',...
                            'Interpreter','latex','Fontsize',fontsize);
            yname = ylabel(ah.YLabel.String,...
                            'Interpreter','latex','Fontsize',fontsize);
            %% Formats legend
            if leg_on
                leg = fh.Children(1);
                set(leg,'Box','off','interpreter','latex');
            end
        % ---- change figure window ---- 
        set(fh,'units',units,'color','white')
        pos_fh = get(fh,'position'); % figure window
        set(fh,'position',[pos_fh(1:2),w1,h1])
        drawnow;
        % ---- change plot axes ----
        if ch_on
            axs = 0.6;
        else
            axs = 0.75; %normalized to width
        end
        set(ah,'units',units,'fontsize',fontsize)
        pos_ah = get(ah,'position'); % plot axes
        set(ah,'position',[pos_ah(1:2),axs*w1,pos_ah(4)])
        % ---- change colorbar axes -----
        if ch_on
            set(ch,'units',units,'fontsize',fontsize)
            pos_ch = get(ch,'position'); % colorbar
            pos_ch(1) = pos_ah(1) + (axs*1.1)*w1;
            set(ch,'position',[pos_ch(1:2),0.5*pos_ch(3),pos_ch(4)])
        end
%         % ---- change ylabel position
%         ynamepos = get(yname,'Position');
%         set(yname,'Position',[(1-2/32)*ynamepos(1), ynamepos(2:3)]);
        % ---- adjust paper size for printing ----
        set(fh,'PaperUnits','centimeters',...
               'PaperPosition',[0,0,w1,h1],...
               'PaperSize',[w1,h1])
        drawnow;
        % Print figures
        if save_on
            disp('Saving...');
                print(fh,fhname,fspec);
                disp(['Plot saved to ',fhname,' as a ',fspec(3:end)]);
        end
        hold off;