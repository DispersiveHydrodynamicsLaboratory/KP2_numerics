function[soli] = zero_mean(soli,Ly,Lx,Ny,x0_odd,tstart)
        soli.y  = linspace(-Ly,Ly-(2*Ly/Ny),Ny);
        soli.dipvec = zeros(size(soli.y));
        for Nyi = 1:length(soli.y)
            soli.dipvec(Nyi) = integral(@(x)gather(soli.ua(x,soli.y(Nyi),tstart)),-Lx,+Lx);
        end
    	soli.x0odd = x0_odd;
        soli.dipw  = 5;
    	soli.dip    = @(x,y)    -1/(pi*soli.dipw)*...
                                interp1(soli.y,soli.dipvec,y).*...
                                sech((x-soli.x0odd)/soli.dipw);
        soli.x0_odd = x0_odd;
        soli.u0    = @(x,y)    soli.ua(x,y,tstart) + soli.dip(x,y);