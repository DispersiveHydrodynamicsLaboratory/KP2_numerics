%% Wrapper function for generic soliton jump
%  with an arbitrary angle
contour_on = 1;
if strcmp(computer,'MACI64')
    updir = '/Users/dhl/Documents/MATLAB/KP_Whitham/';
    slant = '/';
else
    updir = 'H:\';
    slant = '\';
end


%% Parameters, physical space
am = 4;   % Amplitude below jump
ap = 9;   % Amplitude above jump

qm = 1/2; % Slope of propagation below jump
qp = 1; % Slope of propagation above jump

m    =  2;    % Slope where jump in a and q occurs
mrec = -1/m; % Perpendicular to ^

%% Check for acceptable parameters
choice = 0;
% Ensures no shock
if (qm <= qp ) && (qm+qp >= mrec)
    disp('Parameters correspond to Option 1');
    choice = 1;
elseif (qm>qp) && (qm+qp < mrec)
    disp('Parameters correspond to Option 2');
    choice = 2;
else      
    disp(['(qm <= qp ): ',num2str((qm <= qp ))]);
    disp(['qm+qp: ',num2str(qm+qp)]);
    disp(['mrec: ',num2str(mrec)]);
    warning('Incorrect parameters chosen');
end
    
%% Convert to \sqrt{a}
sam = sqrt(am);
sap = sqrt(ap);

%% Compute R1 and R2 on both sides of jump
R1m = sam + qm;
R1p = sap + qp;
R2m = sam - qm;
R2p = sap - qp;

% Ensures 1- and 2- waves are ordered
if ~(R1p^2-2/m*R1p <= R2m^2+2/m*R2m)
    disp(['R1p: ',num2str(R1p),'     R2m: ',num2str(R2m),...
        '     2/m: ',num2str(2/m)]);
    disp(['R2m+R1p:     ',num2str((R2m+R1p    )),...
     '     R2m-R1p+2/m: ',num2str((R2m-R1p+2/m))]);
    if choice == 1
        incfac = (qm+qp-2/m+sap)^2;
        warning(['Warning, wave speeds unphysical.',...
            'Increase am to over ',num2str(incfac),'.']);
    else
        incfac = (-qm-qp-sap+2/m)^2;
        warning(['Warning, wave speeds unphysical.',...
            'Increase am to over ',num2str(incfac),'.']);
    end
end


% Compute solution in terms of Riemann Invariants
[R1,R2,f,g] = KP2jump_angle_no_mean(R1m,R1p,R2m,R2p,m,choice);

% Convert back to physical space
sa = @(x,y,t) 1/2 * (R1(x,y,t) + R2(x,y,t));
a  = @(x,y,t) sa(x,y,t).^2;
q  = @(x,y,t) 1/2 * (R1(x,y,t) - R2(x,y,t));
%% Set up soliton function
u = @(theta,x,y,t,a) a(x,y,t).*(sech(sqrt(a(x,y,t)/12).*theta)).^2;
theta = @(x,y,t,a,q) (x + q(x,y,t).*y - (a(x,y,t)/3+q(x,y,t).^2) * t);


%% Plot Results 
if contour_on
        figname = 'RW_arb_angle';
        fhname = [pwd,slant,figname,'_contour'];
        xmin = -30;   % Set to 0 to begin as close to nozzle as possible
        xmax =  30; % Set to Inf to maximize upper conduit
        ymin = -50;   % Set to 0 to use first time
        ymax =  30; % Set to Inf to use last time 
        t    =   5;
        
        x = xmin:xmax;
        y = ymin:ymax;
        [X,Y] = meshgrid(x,y);
        THETA = theta(X,Y,t,a,q);
        U = u(THETA,X,Y,t,a);

        cmin = 0;
        cred = 0;
        cmax = max(U(:)); % Set to Inf to maximize conduit area

[fh]   = make_contours_for_pub(U,x,y,ymin,ymax,...
                            xmin,xmax,cmin,cred,cmax,...
                            1,fhname,slant);

end