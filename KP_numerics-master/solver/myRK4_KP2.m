function myRK4_KP2( Vhat_init, uasy, dxuasy, dyuasy,...
                    dt, W, Wp, Wpp,...
                    iphi, domain)
%myRK4_KP2 is a modified fourth-order explicit Runge-Kutta timestepping method
% Includes variables and saving methods specific to a KP2 solver
% IMPORTANT: assumes function is not t-dependent
% INPUTS:
%   u_init:      Initial condition
%   v_init:      Initial condition, windowed
%   Vhat_init:   Initial condition, windowed, FFT'd, and int-factored
%   uasy:        Function that matches the unwindowed solution for large y
    %   dxuasy:     its derivative wrt x
    %   dyuasy:     its derivative wrt y
%   dt:          Time stepping increment (approximate)
%   tout:        Times to be output
%   W, Wp, Wpp:  Windowing function and its derivatives (exact)
%   iphi:        Integrating factor exponent
%   domain:      Structure containing the following
    %    x,  y:     real space domain (vectors)
    %    X,  Y:     real space domain (matrices)
    %   dx, dy:     Real space discretization in x and y
    %   Lx, Ly:     Real, unscaled space maxima
    %    k,  l:     wavenumber domain (vectors)
    %   KX, KY:     wavenumber domain (matrices)
% OUTPUTS:
%   (none except to file)

    % Call global variables; initialize RK4 method
    global tout inc dir
    Vold = Vhat_init; % fft'd, unshifted, decayed IC

    % Subsequent time steps
    for jj = 1:length(tout)-1
        disp(['Calculating ',num2str(jj),' out of ',num2str(length(tout)-1)]);
            tmid = linspace(tout(jj),tout(jj+1),ceil((tout(jj+1)-tout(jj))/dt)+1);
            for ii = 2:length(tmid)
                Vnew = RK4(tmid(ii-1), tmid(ii)-tmid(ii-1), Vold, uasy, dxuasy, dyuasy,...
                            W, Wp, Wpp, iphi, domain );
                if sum(isnan(Vnew(:)))>0
                    error(['Not a Number encountered at t=',num2str(tmid(ii))]);
                end
                Vold = Vnew;
            end
        
        %% Save data
        v    = Vnew.*exp(-iphi*tout(jj+1));
        tnow = tout(jj+1);
        u = real(ifft2(v));
          save(strcat(dir,num2str(inc,'%05d')),'u','v','tnow','inc');
          inc = inc +1;
    end


% KP2 RK4 function
function Vhatnew = RK4( t, dt, Vhat, uasy, dxuasy, dyuasy, W, Wp, Wpp, iphi, domain );
% Solves: KP eq. (u_t + uu_x + epsilon^2 u_xxx)_x + lambda u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with integrating factor 
% v = exp[+i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat
%     g  = -1/2*1i*dt*KX;
    % Domain names for ease
    X = domain.X; Y = domain.Y;
    % Exponentials to change Vhat to vhat
    Ezero = exp(-t.*iphi); Ehalf  = exp(-(t+dt/2).*iphi); %Eone = exp(-(t+dt)*iphi);
    vhat = Vhat.*Ezero;
    Va  = G( t     , vhat                ,...
             uasy(X,Y,t)     , dxuasy(X,Y,t)     , dyuasy(X,Y,t)     ,...
             W, Wp, Wpp, iphi, domain );
    Vb  = G( t+dt/2, vhat+Ezero.*Va/2*dt ,...
             uasy(X,Y,t+dt/2), dxuasy(X,Y,t+dt/2), dyuasy(X,Y,t+dt/2),...
             W, Wp, Wpp, iphi, domain );     % 4th-order
    Vc  = G( t+dt/2, vhat+Ehalf.*Vb/2*dt ,...
             uasy(X,Y,t+dt/2), dxuasy(X,Y,t+dt/2), dyuasy(X,Y,t+dt/2),...
             W, Wp, Wpp, iphi, domain );     % Runge-Kutta
    Vd  = G( t+dt  , vhat+Ehalf.*Vc/2*dt ,...
             uasy(X,Y,t+dt)  , dxuasy(X,Y,t+dt)  , dyuasy(X,Y,t+dt)  ,...
             W, Wp, Wpp, iphi, domain );
    Vhatnew = Vhat + dt*(Va + 2*(Vb+Vc) + Vd)/6;
    disp('')
    

function Gv = G( t, v, uasy, dxuasy, dyuasy, W, Wp, Wpp, iphi, domain )
    RHS = (pi/domain.Lx)*(1-W(domain.Y)) .*...
          ifft2(1i*domain.KX.*fft2(W(domain.Y).*uasy.* dxuasy -...
                            ifft2(1i*domain.KX.*fft2(ifft2(v).*uasy)) ))  + ...
          (pi*domain.Lx/domain.Ly) .*...
          ( 2*Wp(domain.Y).*dyuasy + Wpp(domain.Y).*uasy ) ;
    RHShat = fft2(RHS);
	v2hat  = fft2(ifft2(v).^2);
    Gv = exp(iphi*t).*( (1./(1i*domain.KX).*RHShat) -...
                        (1i*domain.KX*pi)./(2*domain.Lx).*v2hat );
	Gv(domain.KX==0) = 0;
    disp('');