function myRK4_KdV( u0, dt, domain)
%myRK4_KdV is a modified fourth-order explicit Runge-Kutta timestepping method
% Includes variables and saving methods specific to a KP2 solver
% IMPORTANT: assumes function is not t-dependent
% INPUTS:
%   u_init:      Initial condition
%   dt:          Time stepping increment (approximate)
%   tout:        Times to be output
%   domain:      Structure containing the following
    %    x:     real space domain (vector)
    %   dx:     Real space discretization in x
    %   Lx:     Real, unscaled space maxima
    %    k:     wavenumber domain (vector)
    
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
                Vnew = RK4(tmid(ii-1), tmid(ii)-tmid(ii-1), Vold, ...
                            u0.ua, u0.uax, u0.uay,...
                            W, iphi, domain );
                if sum(isnan(Vnew(:)))>0
                    error(['Not a Number encountered at t=',num2str(tmid(ii))]);
                end
                % Figure for debugging
%                 if tout(jj+1)>0.5
%                     plot_interim_contours;
%                     disp('');
%                 end
                Vold = Vnew;
            end
        disp('');
        %% Save data
        v    = Vnew.*exp(-iphi*tout(jj+1));
        tnow = gather(tout(jj+1));
        u = gather(real(ifft2(v)) + (1-W.o).*u0.ua(domain.X,domain.Y,tnow));
        v = gather(v);
          save(strcat(dir,num2str(inc,'%05d')),'u','v','tnow','inc');
          inc = inc +1;
    end


% KP2 RK4 function
function Vhatnew = RK4( t, dt, Vhat, uasy, dxuasy, dyuasy, W, iphi, domain );
% Solves: KP eq. (u_t + uu_x + epsilon^2 u_xxx)_x + lambda u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with integrating factor
% v = exp[+i(k^3*epsilon^2-lambda*l^2/k)t]*u_hat
%     g  = -1/2*1i*dt*KX;
    % Domain names for ease
    X = domain.X; Y = domain.Y;
    %% Precompute matrices for speed
    % Exponentials to change Vhat to vhat
    Ezero = exp( t      .*iphi);  Ezeroi = exp(- t      .*iphi);
    Ehalf = exp((t+dt/2).*iphi);  Ehalfi = exp(-(t+dt/2).*iphi);
    Eone  = exp((t+dt)  .*iphi);  Eonei  = exp(-(t+dt  ).*iphi);
    % Function evals of asymptotic solution
        % t = 0
        uazero = uasy(X,Y,t);
        uaxzero = dxuasy(X,Y,t);
        uayzero = dyuasy(X,Y,t);
        uahalf = uasy(X,Y,t+dt/2);
        uaxhalf = dxuasy(X,Y,t+dt/2);
        uayhalf = dyuasy(X,Y,t+dt/2);
        uaone = uasy(X,Y,t+dt);
        uaxone = dxuasy(X,Y,t+dt);
        uayone = dyuasy(X,Y,t+dt);
    Va  = G( Ezero, Ezeroi.* Vhat          ,...
             uazero, uaxzero, uayzero     ,...
             W, domain );
    Vb  = G( Ehalf, Ehalfi.*(Vhat+dt/2*Va) ,...
             uahalf, uaxhalf, uayhalf,...
             W, domain );     % 4th-order
    Vc  = G( Ehalf, Ehalfi.*(Vhat+dt/2*Vb) ,...
             uahalf, uaxhalf, uayhalf,...
             W, domain );     % Runge-Kutta
    Vd  = G( Eone, Eonei  .*(Vhat+dt*Vc) ,...
             uaone  , uaxone  , uayone  ,...
             W, domain );
    Vhatnew = Vhat + dt*(Va + 2*(Vb+Vc) + Vd)/6;
%     plot_interim_contours;
%     disp('');
    

function GV = G( Et, vhat, uasy, dxuasy, dyuasy, W, domain )
    v = ifft2(vhat);
    RHS = (1-W.o) .*...
          ifft2( 1i*domain.KX.*fft2( W.o.*uasy.* dxuasy -...
                 ifft2(1i*domain.KX.*fft2(v.*uasy)) ))  + ...
          ( 2*W.p.*dyuasy + W.pp.*uasy ) ;
    RHShat = fft2(RHS);
	v2hat  = fft2(v.^2);
    GV  = Et.*( 1./(1i*domain.KX).*RHShat -...
                    1i*domain.KX./2.*v2hat );
% 	GV  = Et.*(1i*domain.KX./2.*v2hat );
% 	Gv0 = ( -(1i*domain.KX*pi)./(2*domain.Lx).*v2hat );
    GV(domain.KX==0) = 0;%Gv0(domain.KX==0);
    disp('');
    
    