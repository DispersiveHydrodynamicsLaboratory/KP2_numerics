function[u0] = zero_mean(u0,Ly,Lx,Ny,x0_odd)
    disp('Calculating zero-mean function...');
        u0.dipvec = zeros(size(u0.y));
        for Nyi = 1:length(u0.y)
		q= integral(@(x)gather(u0.us(x,u0.y(Nyi))),-Lx,+Lx,'ArrayValued',1);
            u0.dipvec(Nyi)=	    q;
        end
    	u0.x0odd = x0_odd;
        u0.dipw  = 5;
    	u0.dip    = @(x,y)    -1/(pi*u0.dipw)*...
                                interp1(u0.y,u0.dipvec,y).*...
                                sech((x-u0.x0odd)/u0.dipw);
        u0.x0_odd = x0_odd;
        u0.u0    = @(x,y)    u0.us(x,y) + u0.dip(x,y);
    disp('Function determined, continuing...');