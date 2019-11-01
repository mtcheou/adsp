function sqrnorm = calcSqrNorm(rho,xi,phase,deltat)

if(deltat==0)
    sqrnorm =1.0;
elseif ( (xi~=0.0) && (rho~=0.0) ) 
    sqrnorm = ((cos(2*phase) + 1)*rho*rho - sin(2*phase)*rho*xi + xi*xi)/(4*rho*rho*rho + 4*rho*xi*xi) - ...
              (rho*rho*(cos(2*phase + 2*deltat*xi) + 1) + xi*xi - rho*xi*sin(2*phase + 2*deltat*xi))*exp(-2*deltat*rho)/ ...
              (4*rho*rho*rho + 4*rho*xi*xi);
elseif ( (xi~=0.0) && (rho==0.0) ) 
    sqrnorm = (deltat/2) - (sin(2*phase) - sin(2*phase + 2*deltat*xi))/(4*xi);
elseif ( (xi==0.0) && (rho~=0.0) ) 
    sqrnorm = (cos(phase)*cos(phase))*(1-exp(-2*rho*deltat))/ ...
              (2*rho);
else
    sqrnorm = deltat*(cos(phase)*cos(phase));
end