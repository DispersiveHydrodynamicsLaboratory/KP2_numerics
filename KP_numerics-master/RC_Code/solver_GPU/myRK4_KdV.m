function myRK4_KdV( u0, dt, domain)
%myRK4_KdV is a modified fourth-order explicit Runge-Kutta timestepping method
% Includes variables and saving methods specific to a KdV solver
% IMPORTANT: assumes function is not t-dependent
% INPUTS:
%   u0:      Initial condition
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
    uold = u0; % fft'd, unshifted, decayed IC

    umat = gpuArray(zeros(length(u0),length(tout)));
    umat(:,1) = u0;
    % Subsequent time steps
    for jj = 1:length(tout)-1
        if length(tout)>1000
            if mod(jj,100)==0
                disp(['Calculating ',num2str(jj),' out of ',num2str(length(tout)-1)]);
            end
        else
            disp(['Calculating ',num2str(jj),' out of ',num2str(length(tout)-1)]);
        end     
            tmid = linspace(tout(jj),tout(jj+1),ceil((tout(jj+1)-tout(jj))/dt)+1);
            for ii = 2:length(tmid)
                unew = RK4(tmid(ii-1), tmid(ii)-tmid(ii-1), uold, ...
                            domain );
                if sum(isnan(unew(:)))>0
                    error(['Not a Number encountered at t=',num2str(tmid(ii))]);
                end
                % Figure for debugging
%                 if tout(jj+1)>0.5
%                     plot_interim_contours;
%                     disp('');
%                 end
                uold = unew;
            end
        disp('');
        %% Save data
        umat(:,jj+1) = unew;
%         tnow = gather(tout(jj+1));
%         u = gather(unew);
%           save(strcat(dir,num2str(inc,'%05d')),'u','tnow','inc');
%           inc = inc +1;
    end
    save(strcat(dir,'matrix'),'umat','tout','inc');


% KP2 RK4 function
function unew = RK4( t, dt, u, domain );
% Solves: KdV eq. u_t + uu_x +  u_xxx = 0
    % Domain names for ease
    %x = domain.x; 
    %% Precompute matrices for speed
    
    Va  = G( (u)        , domain );
    Vb  = G( (u+dt/2*Va), domain );     % 4th-order
    Vc  = G( (u+dt/2*Vb), domain );     % Runge-Kutta
    Vd  = G( (u+  dt*Vc), domain );
    unew = u + dt*(Va + 2*(Vb+Vc) + Vd)/6;
%     plot_interim_contours;
%     disp('');
    

function ut = G( u, domain )
% Needs adjusting    
ut = ifft( -1i*domain.kx.*fft(1/2*u.^2) + 1i*(domain.kx.^3).*fft(u)) ;
    