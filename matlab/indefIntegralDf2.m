function f = indefIntegralDf2(rho,xi,phase,t,drho,dphi)

if (rho~=0.0)
    f = ( exp(-2*rho*t) / (4*(rho*rho+xi*xi)) ) * ...
            (xi*sin(2*xi*t+2*phase) - rho* cos(2*xi*t+2*phase) ) + ...
        ( exp(2*(drho-rho)*t)/ (4*(drho-rho)*(drho-rho)+4*xi*xi) ) *  ...
            (xi*sin(2*xi*t+2*phase) + (drho-rho)* cos(2*xi*t+2*phase)) - ...
        ( exp((drho-2*rho)*t)/ ( (drho-2*rho)*(drho-2*rho)+4*xi*xi ) ) * ...
            (xi*sin(2*xi*t+2*phase) + (drho-2*rho)* cos(2*xi*t+2*phase)) + ...
        ( 2*xi*dphi*exp(2*(drho-rho)*t)/ (4*(drho-rho)*(drho-rho)+4*xi*xi) ) * ... 
            (2*(drho-rho)*sin(2*xi*t+2*phase) + 2*xi* cos(2*xi*t+2*phase)) - ...
        ( 2*xi*dphi*exp((drho-2*rho)*t)/ ( (drho-2*rho)*(drho-2*rho)+4*xi*xi ) ) * ...
            ((drho-2*rho)*sin(2*xi*t+2*phase) + 2*xi* cos(2*xi*t+2*phase));
else
    f = 0.0;
end
