function f =  indefIntegralDf1(rho, xi, phase, t, drho, dphi)

if (rho~=0)

    f =  - (exp(-2*rho*t)/(4*rho)) + ( exp(2*(drho-rho)*t)/ (4*(drho-rho)) ) - ...
            cos(dphi)* ( exp((drho-2*rho)*t)/ (drho-2*rho) );
else
    f = t - cos(dphi)*t;
end

