function[x,y,t,U] = RW1wave_evolve_fun(sam,qm,sap,qp,...
                                       xmin,xmax,ymin,ymax,t)

% 1-Wave time evolution
plot_on = 0;

%% Soliton Parameters
if qp ~= qm - sam + sap;
    disp('1-RW conditoins not satisfied');
    qp = qm - sam + sap;
    disp(['Using qp = ',num2str(qp)]);
end


%% Check to make sure characteristic speeds are ok
if  sign(qm)~=sign(qp) ||...
        ~( (0<sam-sap && sam-sap<qm) || (qm<sam-sap && sam-sap<0) )
    disp(['sam - sap: ',num2str(sam-sap)]);
    disp(['qm: ',num2str(qm)]);
    disp(['sign(qm)==sign(qp): ',num2str(sign(qm)==sign(qp))]);
    disp(['(0<sam-sap && sam-sap<qm): ',num2str((0<sam-sap && sam-sap<qm))]);
    disp(['(qm<sam-sap && sam-sap<0): ',num2str((qm<sam-sap && sam-sap<0))]);
    error('Compatibility conditions not qualified');
end

%% Determine if R2 is constant over jump
R2p  = sap - qp;
R2m  = sam - qm;
    if R2p==R2m
        R2 = R2p;
    else
        disp('R2 constant not satisfied');
        return;
    end

%% Work for R1
R1p = sap + qp;
R1m = sam + qm;
l1m = R1m*(2/3*R2 - 1/3*R1m);
l1p = R1p*(2/3*R2 - 1/3*R1p);

R1 = @(x,t)  (R2 + sqrt( R2.^2 - 3*(x/t) )) .* ( x <  l1p*t).*( x >  l1m*t) + ...
             R1m                            .*                ( x <= l1m*t) + ...
             R1p                            .* ( x >= l1p*t);

solap   = @(x,t) ((R1(x,t) + R2)/2).^2;
solqp   = @(x,t)  (R1(x,t)-R2)/2;

%% Set up meshgrid
    x = linspace(xmin,xmax,1000);
    y = linspace(ymin,ymax,1000);
    [X,Y] = meshgrid(x,y);
    xi = x/t;

%% Set up soliton function
u = @(theta,x,t,a) a(x,t).*(sech(sqrt(a(x,t)/12).*theta)).^2;
theta = @(x,y,t,a,q) (x + q(x,t).*y - (a(x,t)/3+q(x,t).^2) * t);
THETA = theta(X,Y,t,solap,solqp);
U = u(THETA,X,t,solap);

if plot_on

    %% Plot Amplitude and angle of propagation; mark characteristic speeds
    figure(2); clf; drawnow;
        subplot(3,1,1);
            plot(x,R1(x,t),'.',[l1m*t, l1p*t],[R1m, R1p],'rx'); xlabel('x'); ylabel('R_1');
            title(['t = ',num2str(t)]);
        subplot(3,1,2); 
            plot(x,solap(x,t),'.'); xlabel('x'); ylabel('a');
        subplot(3,1,3); 
            plot(x,solqp(x,t),'.'); xlabel('x'); ylabel('q');

    %% Plot Full Soliton Solution
    figure(1); clf; drawnow;
    contourf(X,Y,U,'edgecolor','none');
        cmap =load('CoolWarmFloat257.csv');
        colormap(cmap);
        colorbar;
        xlabel('x'); ylabel('y'); drawnow;

end