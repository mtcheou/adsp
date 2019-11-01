function f =  indefIntegralDfb(rho,rhoq, xi, phi, phiq, t)

if ( (xi~=0) && ( (rho+rhoq)~=0) )
    f = - exp(-rho*t - rhoq*t)*(rho*cos(phi + phiq + 2*t*xi) + rhoq*cos(phi + phiq + 2*t*xi) - 2*xi*sin(phi + phiq + 2*t*xi))/ ...
         (((rho*rho) + 2*rho*rhoq + (rhoq*rhoq) + 4*(xi*xi)));
elseif ( (xi~=0) && ( (rho+rhoq)==0) )
    f = sin(phi + phiq + 2*t*xi)/(2*xi);
elseif ( (xi==0) && ( (rho+rhoq)~=0) )
    f = - exp(-rho*t - rhoq*t)*cos(phi + phiq)/(rho + rhoq);
else
    f = t*cos(phi + phiq);
end	