function f = indefIntegralAtomSqrNorm(rho,xi, phase, t)
	
if ( (xi~=0) && (rho~=0) )
    f =- exp(-2*rho*t)*(rho*rho*(cos(2*phase + 2*t*xi) + 1) + (xi*xi) - rho*xi*sin(2*phase + 2*t*xi))/(4*(rho*rho*rho) + 4*rho*(xi*xi));
elseif ( (xi~=0) && (rho==0) )
    f = t/2 + sin(2*phase + 2*t*xi)/(4*xi);
elseif ( (xi==0) && (rho~=0) )
    f = -exp(-2*rho*t)*(cos(phase)*cos(phase))/(2*rho);
else
    f = t*cos(phase)*cos(phase);
end	
