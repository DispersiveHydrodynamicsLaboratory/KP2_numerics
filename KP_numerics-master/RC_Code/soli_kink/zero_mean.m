function[soli] = zero_mean(soli,Ly,Lx,Ny,x0_odd)
        soli.dipvec = zeros(size(soli.y));
	class(soli.y)
        for Nyi = 1:length(soli.y)
		q= integral(@(x)gather(soli.us(x,soli.y(Nyi))),-Lx,+Lx,'ArrayValued',1);
            soli.dipvec(Nyi)=	    q;
        end
    	soli.x0odd = x0_odd;
        soli.dipw  = 5;
    	soli.dip    = @(x,y)    -1/(pi*soli.dipw)*...
                                interp1(soli.y,soli.dipvec,y).*...
                                sech((x-soli.x0odd)/soli.dipw);
        soli.x0_odd = x0_odd;
        soli.u0    = @(x,y)    soli.us(x,y) + soli.dip(x,y);
	class(soli.u0);
	soli.u0(0,0);
