% genexp
function [realAtom normr] = genexp(sb,signalSize)

rho = sb(1);
xi  = sb(2);
phi = sb(3);
a   = sb(4);
b   = sb(5);

normr=0;
realAtom = zeros(1,signalSize);

peak=0;
if (rho>=0)
    u = a;
	for i = 0:signalSize-1
		realAtom(i+1)=0;
		if(i >= u)
			if(xi~=0)
				realAtom(i+1) = exp( -rho * (i-u) ) * cos ( (xi*i) + phi );
			else
				realAtom(i+1) = exp( -rho * (i-u) ) * cos ( phi );
            end
            if ( abs(peak) < abs(realAtom(i+1)) ) 
                peak = realAtom(i+1);
            end
        end
    end			
else
	u = signalSize - 1 - b;
	for i=0: signalSize-1
		realAtom(signalSize-i)=0;
		if( i >= u )
            if(xi~=0)
				realAtom(signalSize-i) = exp( rho * (i-u) ) * cos ( (xi*(signalSize-1-i)) + phi );
			else
				realAtom(signalSize-i) = exp( rho * (i-u) ) * cos ( phi );
            end
			if ( abs(peak) < abs(realAtom(signalSize-i)) ) 
                peak = realAtom(signalSize-i);
            end
        end
    end
end
if (abs(peak) < 1e-10) 
    for i = 0: signalSize-1 
        realAtom(i+1) = 0;
    end
else
    for i=0: a-1
        realAtom(i+1) = 0;
    end
    for i = signalSize - 1:-1: b+1
        realAtom(i+1) = 0;
    end
    if (norm(realAtom)~=0)
        normr = norm(realAtom);
        realAtom = realAtom/norm(realAtom);
    end
end