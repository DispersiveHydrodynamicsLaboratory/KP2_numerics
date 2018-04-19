 am = 1;
sam = sqrt(am);
qm = 0;
 a0 = 0.75;
sa0 = sqrt(a0);
 qh = sa0 - sam;
 ql = sam - sa0;
 
% Plot hodograph variables as functions of (a,q)
avec = 0:0.01:1;
T = 3./(8*avec.^(3/2)).*(1+avec-qm.^2);

figure(1);
    subplot(2,2,1);
        plot(avec,T);
    subplot(2,2,2);
 

