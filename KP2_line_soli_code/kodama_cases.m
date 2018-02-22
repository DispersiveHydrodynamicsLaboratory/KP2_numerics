
% kcase.a = [-3/4 -3/4 1/4 5/4];

%% Values from Kodama's paper
klabel = [{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'},...
          {'g_l'}, {'g_r'}, {'h_l'}, {'h_r'}, {'i_l'}, {'i_r'}, {'j_l'}, {'j_r'}];
Am     = 2*ones(size(klabel));
Ap     = [1/2 49/8 2   2  1/4 6 ...
            2   2   2   2   8/25    8/25    169/32  169/32];
Qm     = [1/2  3/4 1 12/5  0  0 ...
          -2.02 2.02  -1  1  -2/5   2/5     -1/4       1/4];
Qp     = -Qm;

%% Transformation to our version
am = Am*6^(2/3);
ap = Ap*6^(2/3);
qm = Qm*(3/sqrt(2))^(1/3);
qp = Qp*(3/sqrt(2))^(1/3);

%% Check for which case we're in
da = sqrt(am) - sqrt(ap);
dq = qp - qm;

avec = linspace(-15,15,1000);
qvec = linspace(-15/2,15/2,1000);
[AV,QV] = meshgrid(avec,qvec);

vinds = 1:6;
xinds = 7:length(klabel);

fsize = 16;
msize = 12;

fh = figure(1); clf; set(gcf,'Color','White');
    contourf(AV,QV,(-(AV<=QV) - 1/2*(AV>=-QV)));
    set(gca,'fontsize',fsize);
    colormap(gray);
%     colorbar;
    hold on;
    title(['V-case, no x-variation']);
    xlabel(['a_- - a_+']);
    ylabel(['q_+ - q_-']);
        text([-19/2,0,0,19/2],[0,-4,4,0],...
        [{'1S, 2RW'},{'all shocks'},{'all RWs'},{'1RW, 2S'}],...
        'Color','w','fontsize',fsize,'HorizontalAlignment','center');
    text(0,-4,'all shocks','fontsize',fsize,'HorizontalAlignment','center');
    offset = ((-1).^mod(1:length(klabel),2))*2;
    plot(da(vinds),dq(vinds),'r.','MarkerSize',msize);
    text(da(vinds),dq(vinds),klabel(vinds),'Color','r','fontsize',fsize,'HorizontalAlignment','right',...
            'VerticalAlignment','bottom');
    drawnow;
    
print('kodama_v','-dpng');


fh = figure(2); clf; set(gcf,'Color','White');
    contourf(AV,QV,(-(AV<=QV) - 1/2*(AV>=-QV)));
    set(gca,'fontsize',20);
    colormap(gray);
%     colorbar;
    hold on;
    title(['X-case, no x-variation']);
    xlabel(['a_- - a_+']);
    ylabel(['q_+ - q_-']);
        text([-15/2,0,0,15/2],[0,-4,4,0],...
        [{'1S, 2RW'},{'all shocks'},{'all RWs'},{'1RW, 2S'}],...
        'Color','w','fontsize',fsize,'HorizontalAlignment','center');
    text(0,-4,'all shocks','fontsize',fsize,'HorizontalAlignment','center');
    offset = ((-1).^mod(1:length(klabel),2))*2;
    plot(da(xinds),dq(xinds),'r.','MarkerSize',msize);
    text(da(xinds),dq(xinds),klabel(xinds),'Color','r','fontsize',fsize,'HorizontalAlignment','right',...
            'VerticalAlignment','bottom');
    drawnow;
    
print('kodama_x','-dpng');