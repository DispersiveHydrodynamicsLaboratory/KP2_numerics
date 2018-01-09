function[x,y,t,U] = RW2wave_evolve_fun(xmin,xmax,ymin,ymax,t)

% 2-Wave time evolution
plot_on = 0;

%% Soliton Parameters
sam  = 0.5;
qm   = 0.5 ;
sap  = 0.75;
qp   = qm + sam - sap;

%% Check to make sure characteristic speeds are ok
if sign(sap-sam)~=sign(qm) || sign(qm)~=sign(qp) ||...
        ~( (0<sap-sam && sap-sam<qm) || (qm<sap-sam && sam-sap<0) )
    error('Compatibility conditions not qualified');
end
%% Determine if R1 is constant over jump
R1p  = sap + qp;
R1m  = sam + qm;
    if R1p==R1m
        R1 = R1p;
    else
        error('R1 constant not satisfied');
    end
    
%% Work for R2    
R2p = sap - qp;
R2m = sam - qm;
l2m = R2m*( 2/3*R1 - 1/3*R2m );
l2p = R2p*( 2/3*R1 - 1/3*R2p );

R2 = @(x,t)  (R1 - sqrt( R1.^2 - 3*(x/t) )) .* ( x <  l2p*t).*( x >  l2m*t) + ...
             R2m                            .*                ( x <= l2m*t) + ...
             R2p                            .* ( x >= l2p*t);
        
solap   = @(x,t) ((R2(x,t) + R1)/2).^2;
solqp   = @(x,t)  (R1 - R2(x,t))/2;

%% Setup meshgrid
    x = linspace(xmin,xmax,1000);
    y = linspace(ymin,ymax,1000);
    [X,Y] = meshgrid(x,y);
    xi = x/t;

%% Set up soliton function
u = @(theta,x,t,a) a(x,t).*(sech(sqrt(a(x,t)/12).*theta)).^2;
theta = @(x,y,t,a,q) x + q(x,t).*y - (a(x,t)/3+q(x,t).^2) * t;
THETA = theta(X,Y,t,solap,solqp);
U = u(THETA,X,t,solap);


if plot_on
    %% Plot Amplitude and angle of propagation
    figure(2); clf; drawnow;
        subplot(3,1,1);
            plot(x,R2(x,t),'.',[l2m*t, l2p*t],[R2m, R2p],'rx'); xlabel('x'); ylabel('R_2');
            title(['t = ',num2str(t)]);
        subplot(3,1,2); 
            plot(x,solap(x,t),'.'); xlabel('x'); ylabel('a');
        subplot(3,1,3); 
            plot(x,solqp(x,t),'.'); xlabel('x'); ylabel('q');

    %% Plot Full Soliton Solution
    figure(1); clf; drawnow;
    contourf(X,Y,real(U),'edgecolor','none');
        cmap =load('CoolWarmFloat257.csv');
        colormap(cmap);
        colorbar;
        xlabel('x'); ylabel('y'); drawnow;
end