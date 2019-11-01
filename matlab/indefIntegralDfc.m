function f =  indefIntegralDfc(rho,rhoq, phi, phiq, t)


if  ((rho+rhoq)~=0)
    f = - exp(-rho*t - rhoq*t)*cos(phi - phiq)/(rho + rhoq);
else
    f = t*cos(phi - phiq);
end	