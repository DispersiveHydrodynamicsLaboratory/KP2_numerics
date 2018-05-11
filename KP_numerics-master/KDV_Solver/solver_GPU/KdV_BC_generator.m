function[u0] = KdV_BC_generator( u0, t, Lx, Nx, Nt, maindir )
%% MM: FUTURE EDIT: put hi and lo as substructures
disp('Determining KdV boundary conditions...');

%% First, find (or run) KdV version of BC conditions
  KdV_dir = [maindir,filesep,'KP',filesep,'KdV',filesep];
  if ~exist(KdV_dir,'dir')
      mkdir(KdV_dir);
  end
   hi_dir = [KdV_dir,num2str(u0.hi,'%05.2f'),filesep];%'.mat'];
   u0.hi_dir = hi_dir;
   lo_dir = [KdV_dir,num2str(u0.lo,'%05.2f'),filesep];%'.mat']
   u0.lo_dir = lo_dir;
  if exist(hi_dir)~=7
        mkdir(hi_dir);
        try
            save([hi_dir,'parameters.mat'],'t','Lx','Nx','u0');
            KdV_solver_periodic(  t, Lx, Nx, Nt, u0.ubhi, hi_dir );
        catch ehi
            disp(ehi.message);
            rmdir(hi_dir,'s');
        end
  end
  if exist(lo_dir)~=7
        mkdir(lo_dir)
        try
            save([lo_dir,'parameters.mat'],'t','Lx','Nx','u0');
            KdV_solver_periodic(  t, Lx, Nx, Nt, u0.ublo, lo_dir );
        catch elo
            disp(elo.message);
            rmdir(lo_dir);
        end
  end
  
  disp('Boundary Conditions found.');
% %% Then, incorporate KdV results into KP BCs
%     % Load KdV results as matrices in (x,t)
%     load([hi_dir,'matrix.mat'],'umat','tout');
%     u0.uasy.u.dip.hi = umat;
%     load([lo_dir,'matrix.mat'],'umat');
%     u0.uasy.u.dip.lo = umat;
%     % Corresponding x derivatives (y-derivs are 0)
%     u0.dx = (2*Lx/Nx);
%     u0.uasy.ux.dip.hi = ( [u0.uasy.u.dip.hi(end,:);u0.uasy.u.dip.hi] -...
%                         [u0.uasy.u.dip.hi(2:end,:); u0.uasy.u.dip.hi(1:2,:)] )./u0.dx;
%     u0.uasy.ux.dip.lo = ( [u0.uasy.u.dip.lo(end,:);u0.uasy.u.dip.lo] -...
%                         [u0.uasy.u.dip.lo(2:end,:); u0.uasy.u.dip.lo(1:2,:)] )./u0.dx;
%     
%     u0.uasy.u.dip.hif  = @(t) interp2(tout,(Lx/(2*Nx))*(-Nx/2:Nx/2-1)',u0.uasy.u.dip.hi,...
%                     t,(Lx/(2*Nx))*(-Nx/2:Nx/2-1)');
%     u0.uasy.ux.dip.hif = @(t) interp2(tout,(Lx/(2*Nx))*(-Nx/2:Nx/2-1)',u0.uasy.ux.dip.hi,...
%                     t,(Lx/(2*Nx))*(-Nx/2:Nx/2-1)');
%     u0.uasy.u.dip.lof  = @(t) interp2(tout,(Lx/(2*Nx))*(-Nx/2:Nx/2-1)',u0.uasy.u.dip.lo,...
%                     t,(Lx/(2*Nx))*(-Nx/2:Nx/2-1)');
%     u0.uasy.ux.dip.lof = @(t) interp2(tout,(Lx/(2*Nx))*(-Nx/2:Nx/2-1)',u0.uasy.ux.dip.lo,...
%                     t,(Lx/(2*Nx))*(-Nx/2:Nx/2-1)');
%     u0.ua  = @(x,y,t)  u0.uasy.u.exa(x,y,t) +  u0.uasy.u.dip.hif(t).*(x>0) +  u0.uasy.u.dip.lof(t).*(x<=0);
% 	u0.uax = @(x,y,t)  u0.uasy.ux.exa(x,y,t) + u0.uasy.ux.dip.hif(t).*(x>0) + u0.uasy.ux.dip.lof(t).*(x<=0);
%     u0.uay = @(x,y,t)  u0.uasy.uy.exa(x,y,t);
%     
%     disp('');
    