%% Trial FFT solver to check scaling, wavenumbers.

Lx = 50;
Nx = 2^9;

x = (2*Lx/Nx)*(-Nx/2:Nx/2-1)';
k = (pi/Lx)*[0:Nx/2-1 0 -Nx/2+1:-1]';

y = 3*sin(x) + 1/2*sin(3*x) + cos(5*x);
yhat = fft(y);

figure(1);
    subplot(2,1,1);
        plot(x,y);
	subplot(2,1,2);
        plot(k,real(yhat),'b.',...
             k,imag(yhat),'r.');