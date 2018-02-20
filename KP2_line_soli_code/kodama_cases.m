klabel = [{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'}];
% kcase.a = [-3/4 -3/4 1/4 5/4];

%% Values from Kodama's paper
Am = 2*ones(size(klabel));
Ap = [1/2 49/8 2   2  1/4 6 ];
Qm = [1/2  3/4 1 12/5  0  0];
Qp = -Qm;

%% Transformation to our version
am = Am*6^(2/3);
ap = Ap*6^(2/3);
qm = Qm*(3/sqrt(2))^(1/3);
qp = Qp*(3/sqrt(2))^(1/3);

%% Check for which case we're in
da = am - ap;
dq = qp - qm;

avec = linspace(-15,15,1000);
qvec = linspace(-15,15,1000);
[AV,QV] = meshgrid(avec,qvec);

fh = figure(1); clf; set(gcf,'Color','White');
    contourf(AV,QV,((AV<=QV) + 1/2*(AV>=-QV)));
    set(gca,'fontsize',20);
    colormap(gray);
%     colorbar;
    hold on;
    text(da,dq,klabel,'Color','white','fontsize',20);
    drawnow;
    
print('kodama_cases','-dpng');
