% genexpsample
function realAtomSample = genexpsample(sb,signalSize,t)

rho = sb(1);
xi  = sb(2);
phi = sb(3);
a   = sb(4);
b   = sb(5);

i=t;
if (rho>=0)
    u = a;
%    if(i >= u)
        if(xi~=0)
            realAtomSample = exp( -rho * (i-u) ) * cos ( (xi*i) + phi );
        else
            realAtomSample = exp( -rho * (i-u) ) * cos ( phi );
        end
%    end    			
else
	u = signalSize-1-b;
%    if( i <= u )
        if(xi~=0)
            realAtomSample = exp( rho * (signalSize-1-i-u) ) * cos ( (xi*(i)) + phi );
        else
            realAtomSample = exp( rho * (signalSize-1-i-u) ) * cos ( phi );
        end
%    end
end