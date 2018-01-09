function myRK4_KP2( iphi, KX,KY, W, dt, tout, yinit, data_dir)
%myRK4_KP2 is a modified fourth-order explicit Runge-Kutta timestepping method
% Includes variables and saving methods specific to a KP2 solver
% IMPORTANT: assumes function is not t-dependent
% INPUTS:
%   odefun: function to be evaluated odefun = odefun(t,y,...)
%   dt:     time stepping increment (approximate)
%   tout:   times to be output
%   (kx,ky):      wavenumber domain (matrices)
%   yinit:  initial condition, FFT
% OUTPUTS:
%   (none except to file)
if tout(1)~=0
    tout = [0 tout];
end

% Prepare for looping over time steps
% Save initial data
	v    = yinit;
    tnow = 0; inc = 0;
    u = real(ifft2(v));
      save(strcat(data_dir,num2str(1,'%05d')),'u','v','tnow','inc');
      inc = inc + 1;

% Subsequent time steps
for jj = 2:length(tout)
	disp(['Calculating ',num2str(jj),' out of ',num2str(length(tout))]);
        tmid = linspace(tout(jj-1),tout(jj),ceil((tout(jj)-tout(jj-1))/dt));
        yold = v;
        for ii = 2:length(tmid)
            ynew = RK4(tmid(ii)-tmid(ii-1), yold, iphi, KX);
            if sum(isnan(ynew))>0
                error(['Not a Number encountered at t=',num2str(tmid(ii))]);
            end
        end
        hold off;
	% Save data
          v1    = ynew;
          tnow = tout(jj); inc = jj;
      u1 = real(ifft2(v1));
      u = u.*W; % Apply window function
%       % Figure for debugging window fcn
%       figure(5); clf;
%         subplot(2,2,1);
%             contourf(u1,100,'edgecolor','none'); colorbar;
%             title('Original');
%         subplot(2,2,2);
%             contourf(u,100,'edgecolor','none'); colorbar;
%             title('Windowed');
%         subplot(2,2,[3 4]);
%             contourf(u1-u,100,'edgecolor','none'); colorbar;
%             title('u1-u');
%             drawnow;
%             input('R');
      v = fft2(u);
      save(strcat(data_dir,num2str(jj,'%05d')),'u','v','tnow','inc');
      inc = inc +1;
end

% if plot_on
%     figure(5);
%     for jj = 1:length(t)
%         subplot(length(t),1,jj);
%             fftplot(x,y(jj,:));
%             hold on;
%     end
%     drawnow;
% end

% KP2 RK4 function
function vt = RK4( dt, v, iphi, KX)
% Solves: KP eq. (u_t + uu_x + epsilon^2 u_xxx)_x + lambda u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with integrating factor 
% v = exp[+i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat
%% TODO: CHANGE COEFF ON uu_x to 1
    g  = -1/2*1i*dt*KX;
    E  = exp(dt/2*iphi); E2 = E.^2;
    a  = g.*fft2(real( ifft2(     v    ) ).^2);
    b  = g.*fft2(real( ifft2(E.*(v+a/2)) ).^2);     % 4th-order
    c  = g.*fft2(real( ifft2(E.*v + b/2) ).^2);     % Runge-Kutta
    d  = g.*fft2(real( ifft2(E2.*v+E.*c) ).^2);
    vt = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;

