
if 0 %Plot transformed ICs
    figure(3); clf;
        subplot(2,3,1);
            contourf(domain.X,domain.Y,u_init,100,'edgecolor','none');
            colorbar;
        subplot(2,3,2);
            contourf(domain.X,domain.Y,ifft2(fft2(v_init)),100,'edgecolor','none');
            colorbar;
        subplot(2,3,3);
            contourf(domain.X,domain.Y,ifft2(Vhat_init),100,'edgecolor','none');
            colorbar;
%         subplot(2,3,4);
%             contourf(domain.KX,domain.KY,abs(Vc),100,'edgecolor','none')
%             colorbar;
%         subplot(2,3,5);
%             contourf(domain.KX,domain.KY,abs(Vd),100,'edgecolor','none')
%             colorbar;
end

if 0 % Plot RK4 IC, beginning and end
    figure(2); clf;
        subplot(2,2,1);
            contourf(domain.X,domain.Y,ifft2(Vhat_init),100,'edgecolor','none');
            colorbar;
        subplot(2,2,2);
            contourf(domain.X,domain.Y,abs(ifft2(Vold)),100,'edgecolor','none');
            colorbar;
        subplot(2,2,3);
            contourf(domain.X,domain.Y,abs(ifft2(Vnew)),100,'edgecolor','none')
            colorbar;
end

if 0 % Plot RK4 intermediates
    figure(3); clf;
        subplot(2,3,1);
            contourf(domain.KX,domain.KY,abs(vhat),100,'edgecolor','none');
            colorbar;
        subplot(2,3,2);
            contourf(domain.KX,domain.KY,abs(Va),100,'edgecolor','none');
            colorbar;
        subplot(2,3,3);
            contourf(domain.KX,domain.KY,abs(Vb),100,'edgecolor','none');
            colorbar;
        subplot(2,3,4);
            contourf(domain.KX,domain.KY,abs(Vc),100,'edgecolor','none')
            colorbar;
        subplot(2,3,5);
            contourf(domain.KX,domain.KY,abs(Vd),100,'edgecolor','none')
            colorbar;
end

if 0 % Troubleshooting RK4 intermediates
        figure(3); clf;
        subplot(2,3,1);
            contourf(domain.KX,domain.KY,abs(vhat),100,'edgecolor','none');
            colorbar;
        subplot(2,3,2);
            contourf(domain.KX,domain.KY,abs(Va),100,'edgecolor','none');
            colorbar;
        subplot(2,3,3);
            contourf(domain.KX,domain.KY,abs(vhat+Ezero.*Va/2*dt),100,'edgecolor','none');
            colorbar;
%         subplot(2,3,4);
%             contourf(domain.KX,domain.KY,abs(Vb),100,'edgecolor','none')
%             colorbar;
end

if 0 % Plot derivative intermediates
    figure(3); clf;
        subplot(2,3,1);
            contourf(domain.X,domain.Y,ifft2(v),100,'edgecolor','none');
            colorbar;
        subplot(2,3,2);
            contourf(domain.X,domain.Y,RHS,100,'edgecolor','none');
            colorbar;
        subplot(2,3,3);
            contourf(domain.KX,domain.KY,abs(RHShat),100,'edgecolor','none');
            colorbar;
        subplot(2,3,4);
            contourf(domain.KX,domain.KY,abs(v2hat),100,'edgecolor','none')
            colorbar;
        subplot(2,3,5);
            contourf(domain.KX,domain.KY,real(Gv),100,'edgecolor','none')
            colorbar;
        subplot(2,3,6);
            contourf(domain.KX,domain.KY,imag(Gv),100,'edgecolor','none')
            colorbar;
end

if 1 % IC and other plots from loading data directory
    numouts = [1 2 3 4];
	load(strcat(data_dir,num2str(0,'%05d'),'.mat'),'u_init');
    load([data_dir,'parameters.mat'],'t','Nx','Ny','Lx','Ly'); %,'f');
        x = (2*Lx/Nx)*(-Nx/2:Nx/2-1)';
        y = (2*Ly/Ny)*(-Ny/2:Ny/2-1)';
    figure(3); clf;
        subplot(2,3,1);
            contourf(x,y,u_init,100,'edgecolor','none');
            colorbar;
    for ii = 1:length(numouts)
        load(strcat(data_dir,num2str(numouts(ii),'%05d'),'.mat'),'u');
        subplot(2,3,ii+1);
            contourf(x,y,u,100,'edgecolor','none');
            colorbar;
    end 
    drawnow;
end
