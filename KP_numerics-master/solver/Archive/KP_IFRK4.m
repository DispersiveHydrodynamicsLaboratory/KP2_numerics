% Solve KP eq. (u_t + 6uu_x + epsilon^2 u_xxx)_x + lambda u_yy = 0
% on [-Lx,Lx] & [-Ly,Ly] by FFT in space with integrating factor 
% v = exp[-i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat and RK4 in time

clear all;
close all;

% Set up grid and the initial data;

Nx = 2^9; Ny = 2^5; dt = 1.e-3;
Lx = 40; Ly = 10;
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

% Choose IC

%% Set up soliton function
a =@(x,y) 1; q =@(x,y) 0;
u = @(theta,x,y,a) a(x,y).*(sech(sqrt(a(x,y)/12).*theta)).^2;
theta = @(x,y,a,q,x0,y0) ( (x-x0) + q(x,y).*(y-y0) );
IC = @(X,Y) u(theta(X,Y,a,q,10,0),X,Y,a);
u0 = IC(X,Y);

% 0.1 KdV soliton
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

u = u0;
v = fft2(u);

% Solve PDE and plot results:

tmax = 1; nmax = round(tmax/dt); nplt = floor((nmax/10));
for n = 1:nmax
    
    t = n*dt;
    fprintf('t = %.4f\n',t)
    g = -3*1i*dt*KX;
    E = exp(dt*ik3/2); E2 = E.^2;
    a = g.*fft2(real( ifft2(     v    ) ).^2);
    b = g.*fft2(real( ifft2(E.*(v+a/2)) ).^2);     % 4th-order
    c = g.*fft2(real( ifft2(E.*v + b/2) ).^2);     % Runge-Kutta
    d = g.*fft2(real( ifft2(E2.*v+E.*c) ).^2);
    v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;

    vn = norm(v,inf);
    fprintf('norm = %.4f\n',vn)
    
    if mod(n,nplt) == 0
        u = real(ifft2(v));
        surf(X,Y,u);
        shading interp 
        axis([-Lx,Lx,-Ly,Ly,-1,1]);
        drawnow; pause(0.25);
    end

end

% Error for the one-soliton solution of the KdV
%   u = real(ifft2(v));
%   uexact = 2*C^2*sech(C*X-4*C^3*t).^2;
%   error = norm(u-uexact,inf);


u = real(ifft2(v));
surf(X,Y,u);
shading interp 
axis([-Lx,Lx,-Ly,Ly,-1,2]);






