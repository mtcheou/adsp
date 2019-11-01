function [df kappa2 kappaq2] = calcAtomDistFunction(sb,sbQ,sigSize)
 
innerProd = sb(1);
rho = sb(2);
xi = sb(3);
phase = sb(4);
a = sb(5);
b = sb(6);

rhoq = sbQ(2);
phaseq = sbQ(4);

if (a~=b)

    if ((rho<0.0))
        rho = -rho;
        phase = - phase;
        rhoq = -rhoq;
        phaseq = -phaseq;
    %     if (phase<0) 
    %         phase =phase+2*pi;
    %     end
        a = sigSize-1 - sb(6);
        b = sigSize-1 - sb(5);
    end

%     sqrnorm3 =   indefIntegralAtomSqrNorm(rho,xi,phase,b-a) - ...
%                 indefIntegralAtomSqrNorm(rho,xi,phase,0);
%             
%     sqrnorm2 =   indefIntegralAtomSqrNorm(rho,xi,phase,b+1) - ...
%                 indefIntegralAtomSqrNorm(rho,xi,phase,a);
%             
%     sqrnorm = calcSqrNorm(rho,xi,phase,(b-a));
    
    kappa2 =  indefIntegralAtomSqrNorm(rho,xi,phase,b-a) - ...
              indefIntegralAtomSqrNorm(rho,xi,phase,0);
    kappaq2 =  indefIntegralAtomSqrNorm(rhoq,xi,phaseq,b-a) - ...
               indefIntegralAtomSqrNorm(rhoq,xi,phaseq,0);

    df = (1/kappa2)  * ( indefIntegralDfa(rho, xi, phase, (b-a)) - ...
                         indefIntegralDfa(rho, xi, phase,  0   )) + ...
         (1/kappaq2) *(  indefIntegralDfa(rhoq, xi, phaseq, (b-a)) - ...
                         indefIntegralDfa(rhoq, xi, phaseq,  0 )) - ...
         (1/(sqrt(kappaq2)* sqrt(kappa2)))* ... 
                        (   indefIntegralDfb(rho,rhoq, xi, phase, phaseq,(b-a)) - ...
                            indefIntegralDfb(rho,rhoq, xi, phase, phaseq, 0   )) - ...
         (1/(sqrt(kappaq2)* sqrt(kappa2)))* ... 
                        (   indefIntegralDfc(rho,rhoq, phase, phaseq,(b-a)) - ...
                            indefIntegralDfc(rho,rhoq, phase, phaseq, 0   ));
    
else
    if ( ( (cos(phase) > 0.0) && ( cos(phaseq) >0.0) ) ||  ... 
         ( (cos(phase) < 0.0) && ( cos(phaseq) <0.0) ) )
        df = 0;
    elseif ( ( ( cos(phase) == 0.0) && ( cos(phaseq) ~= 0.0) ) ||  ... 
             ( ( cos(phase) ~= 0.0) && ( cos(phaseq) == 0.0) ) )
        df = 1;
    else
        df = 4; 
    end
    
    kappa2 = 1;
    kappaq2 = 1;
    
end

