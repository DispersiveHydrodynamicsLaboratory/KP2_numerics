function vt = compute_deriv( t, v, iphi, KX)
% Solves: KP eq. (u_t + 6uu_x + epsilon^2 u_xxx)_x + lambda u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with integrating factor 
% v = exp[+i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat and Matlab's ode45.m in time
% vt = exp(iphi*t) * 3*i*KX * fft2( ( ifft2( exp(-iphi*t)*v ) )^2 );

 vt = exp(iphi*t) .* 3.*1i.*KX .* fft2( ( ifft2( exp(-iphi*t).*v ) ).^2 );