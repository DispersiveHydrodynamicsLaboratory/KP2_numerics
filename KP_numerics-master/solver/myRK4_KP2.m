function myRK4_KP2( Vhat_init, uasy, dt, tout, W, Wp, Wpp,...
                    iphi, D1x, D1y, KX)
%myRK4_KP2 is a modified fourth-order explicit Runge-Kutta timestepping method
% Includes variables and saving methods specific to a KP2 solver
% IMPORTANT: assumes function is not t-dependent
% INPUTS:
%   u_init:     Initial condition
%   v_init:     Initial condition, windowed
%   Vhat_init:  Initial condition, windowed, FFT'd, and int-factored
%   uasy:       Function that matches the unwindowed solution for large y
%   dt:         Time stepping increment (approximate)
%   tout:       Times to be output
%   W, Wp, Wpp: Windowing function and its derivatives (exact)
%   iphi:       Integrating factor exponent
%   D1, D2:     2nd order Finite Difference Matrices for 1st and 2nd derivs
%   KX, Ky:     wavenumber domain (matrices)
% OUTPUTS:
%   (none except to file)

% Initialize solver
inc = 1;
Vold = Vhat_init;


% Subsequent time steps
for jj = 2:length(tout)
	disp(['Calculating ',num2str(jj),' out of ',num2str(length(tout))]);
        tmid = linspace(tout(jj-1),tout(jj),ceil((tout(jj)-tout(jj-1))/dt));
        for ii = 2:length(tmid)
            Vnew = RK4(tmid, tmid(ii)-tmid(ii-1), Vold, uasy,...
                        W, Wp, Wpp, iphi, D1x, D1y, KX );
            if sum(isnan(Vnew))>0
                error(['Not a Number encountered at t=',num2str(tmid(ii))]);
            end
            Vold = Vnew;
        end
	%% MM: BELOW NEEDS EDITS
	%% Save data
	v1    = Vnew;
	tnow = tout(jj); inc = jj;
	u = real(ifft2(v1));
	u = u.*W; % Apply window function
      %% Save v
      v = fft2(u);
      save(strcat(data_dir,num2str(jj,'%05d')),'u','v','tnow','inc');
      inc = inc +1;
end
end

% KP2 RK4 function
function Vhatnew = RK4( t, dt, Vhat, uasy, W, Wp, Wpp, iphi, D1x, D1y, KX );
% Solves: KP eq. (u_t + uu_x + epsilon^2 u_xxx)_x + lambda u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with integrating factor 
% v = exp[+i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat
%     g  = -1/2*1i*dt*KX;
    % Exponentials to change Vhat to vhat
    Ezero = exp(-t*iphi); Ehalf  = exp(-(t+dt/2)*iphi); %Eone = exp(-(t+dt)*iphi);
    vhat = Vhat*Ezero;
    Va  = G( t     , vhat               , uasy(X,Y,t)     , W, Wp, Wpp, iphi, D1x, D1y, KX );
    Vb  = G( t+dt/2, vhat+Ezero*Va/2*dt , uasy(X,Y,t+dt/2), W, Wp, Wpp, iphi, D1x, D1y, KX );     % 4th-order
    Vc  = G( t+dt/2, vhat+Ehalf*Vb/2*dt , uasy(X,Y,t+dt/2), W, Wp, Wpp, iphi, D1x, D1y, KX );     % Runge-Kutta
    Vd  = G( t+dt  , vhat+Ehalf*Vc/2*dt , uasy(X,Y,t+dt)  , W, Wp, Wpp, iphi, D1x, D1y, KX );
    Vhatnew = Vhat + dt*(Va + 2*(Vb+Vc) + Vd)/6;

end    

function Gv = G( t, v, uasy, W, Wp, Wpp, iphi, D1x, D1y, KX )
    RHS = (pi/Lx)*(1-W(Y)) .* (D1x*( W(Y).*uasy.*(D1x*uasy) - ...
                               D1x*(v.*uasy) )) + ...
          (pi*Lx/Ly) * ( 2*Wp(Y).*(D1y*uasy + Wpp(Y)*uasy ) );
    RHShat = fft2(RHS);
	v2hat  = fft2(v.^2);
    Gv = exp(iphi*t)*( 1./(1i*KX).*RHShat - (1i*KX*pi)./(2*Lx).*v2hat );
end
