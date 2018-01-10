% Solve KP eq. (u_t + 6uu_x + epsilon^2 u_xxx)_x + lambda u_yy = 0
% on [-Lx,Lx] & [-Ly,Ly] by FFT in space with integrating factor 
% v = exp[-i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat and RK4 in time
rip;

% Set up grid and the initial data;
Nx = 2^9; Ny = 2^5; dt = 1.e-3;
Lx = 100; Ly = 50;
% lambda = -1;    % KP I
lambda = 1;     % KP II
o = 1.e-16; epsilon = 1;

x = (2*Lx/Nx)*(-Nx/2:Nx/2-1)';
y = (2*Ly/Ny)*(-Ny/2:Ny/2-1)';
[X,Y] = meshgrid(x,y);
k = (pi/Lx)*[0:Nx/2-1 0 -Nx/2+1:-1]';
l = (pi/Ly)*[0:Ny/2-1 0 -Ny/2+1:-1]';
[KX,KY] = meshgrid(k,l);
ik3 = 1i*(epsilon^2*KX.^3-lambda*KY.^2./(KX+1i*lambda*o));

% Solve PDE and plot results:

tmax = 5; nmax = round(tmax/dt); 

% If plotting, number of times to plot
nplt = 6; 
pltout = round(linspace(1,nmax,nplt));
t = pltout;
figure(2); clf;


% Choose IC

%% Set up soliton function
a =@(x,y) 3; q =@(x,y) 1/4;
u = @(theta,x,y,a) a(x,y).*(sech(sqrt(a(x,y)/12).*theta)).^2;
theta = @(x,y,a,q,x0,y0) ( (x-x0) + q(x,y).*(y-y0) );
IC = @(X,Y) u(theta(X,Y,a,q,25,0),X,Y,a);
% u0 = IC(X,Y)-IC(-X,Y);
ic_type = 'soliton';

%% Factor zeroing IC as y->ymax
winspan = 2.5; % width of decay window
w = @(Y) repmat(tukeywin(size(Y,1),winspan/Ly),[1,size(Y,2)]);

u0 = (IC( X,Y)-IC(X+Lx/2,Y)).*w(Y);
% u0 = IC(X,Y).*w(Y);

%% 0.1 KdV soliton
% C = 1; u0 = 2*C^2*sech(C*X).^2;   % epsilon = 1;

% 0.2 cnidal wave
% m = 0.5;  % r_1 = 0; r_2 = m; r_3 = 1;
% [SN,CN,DN] = ellipj(X,m);   % epsilon = 1;
% u0 = -m+1+2*m*CN.^2;

% 1. IC1 (Klein)
% u0 = -sech(sqrt(X.^2+Y.^2)).^2;
% u0_hat = fft2(u0);
% w_hat = 1i*KX.*u0_hat;
% u0 = real(ifft2(w_hat));

% 2. Parabolic front IC
% mu = -1; 
% u0 = 0.5*(mu*tanh(10*(X+0.01*Y.^2/2))-mu*tanh(10*(X+30+0.01*Y.^2/2)));  

% 3. Cos front IC
% mu = -1;
% u0 = 0.5*(mu*tanh(10*(X-10-cos(pi*Y/Ly)))-mu*tanh(10*(X+10-0.01*cos(pi*Y/Ly))));  

if strcmp(computer,'MACI64')
    maindir = '/Users/dhl/Documents/MATLAB';
    slant = '/';
else
    maindir = 'H:';
    slant = '\';
end

%% Create directory run will be saved to
data_dir = [maindir,slant,'Numerics',slant,'KP',slant,...
            '_tmax_',   num2str(round(tmax)),...
            '_Lx_',   num2str(round(Lx)),...
            '_Nx_',     num2str(Nx),...
            '_Ly_',   num2str(round(Ly)),...
            '_Ny_',     num2str(Ny),...
            '_epsilon_',num2str(epsilon),...
            '_lambda_', num2str(lambda),...
            '_bndry_condns_','periodic',...
            '_init_condns_',ic_type,...
            slant];
% Create the data directory if necessary
if ~exist(data_dir,'dir')
    mkdir(data_dir);
else
    disp(['Warning, directory ',data_dir]);
    disp('already exists, possibly overwriting data');
end

savefile = sprintf('%sparameters.mat',data_dir);
save(savefile);

%% Set up and save IC
u = u0;
v = fft2(u);
inc = 0;
tnow = 0;
save(strcat(data_dir,num2str(inc,'%05d')),'u','v','tnow','inc');
inc = inc + 1;

tic
for n = 1:nmax
    
    t = n*dt;
    
    g = -1/2*1i*dt*KX;
    E = exp(dt*ik3/2); E2 = E.^2;
    a = g.*fft2(real( ifft2(     v    ) ).^2);
    b = g.*fft2(real( ifft2(E.*(v+a/2)) ).^2);     % 4th-order
    c = g.*fft2(real( ifft2(E.*v + b/2) ).^2);     % Runge-Kutta
    d = g.*fft2(real( ifft2(E2.*v+E.*c) ).^2);
    v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
    
    if ismember(n,pltout)
        tnow = t; 
        fprintf('t = %.4f\n',t)
            vn = norm(v,inf);
            fprintf('norm = %.4f\n',vn)
        figure(2);
            h((pltout==n))=subplot(nplt,1,find(pltout==n));
                u = real(ifft2(v));
                contourf(X,Y,u,100,'edgecolor','none');
                colormap(load('CoolWarmFloat257.csv'));
                axis([-Lx,Lx,-Ly,Ly]);
                caxis([-1.03 1.0383]);
                title(['t = :',num2str(t,'%.4f')]);
        drawnow; pause(0.25);
        save(strcat(data_dir,num2str(inc,'%05d')),'u','v','tnow','inc');
        inc = inc + 1;
        toc;
        %% Reapply window function to u
        u = u.*w(Y);
        v = fft2(u);
    end

end


% Error for the one-soliton solution of the KdV
%   u = real(ifft2(v));
%   uexact = 2*C^2*sech(C*X-4*C^3*t).^2;
%   error = norm(u-uexact,inf);


% u = real(ifft2(v));
% surf(X,Y,u);
% shading interp 
% axis([-Lx,Lx,-Ly,Ly,-1,2]);






