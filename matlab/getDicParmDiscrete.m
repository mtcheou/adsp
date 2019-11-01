function dicParm = getDicParmDiscrete(options,N)
%
% options - 'std', 'noFreq', 'noTimeShift'
%
%

%%%%%%%%%%%%%
L=log2(N);

a=2;
l=1;
P_FIM = zeros(1,L+1); 
K_FIM = zeros(1,L+1);
numAtom=0;

if strcmp(options,'std')
    for j=0:L
        P_FIM(j+1) = N*(a^(l-1-j))-1;
        K_FIM(j+1) = ceil((a^(l-1+j))/2) - 1;
        for p=0: P_FIM(j+1)
            for k=0: K_FIM(j+1)
                numAtom=numAtom+1;
            end
        end
    end
    dicParm = zeros(numAtom,3);
    indAtom=0;
    for j=0:L
        s = a^j;
        %s=2^j;
        for p=0: P_FIM(j+1)
            u = p*(a^(j+1-l));
            %u=p*(2^j);
            for k=0: K_FIM(j+1)
                xi = k*pi*a^(2-j-l);
                %xi = k*pi*2^(1-j);
                indAtom=indAtom+1;
                dicParm(indAtom,:) = [s u xi];
            end
        end
    end
elseif strcmp(options,'freqUnif')
    for j=0:L
        P_FIM(j+1) = N*(a^(l-1-j))-1;
        K_FIM(j+1) = (N/2) - 1;
        for p=0: P_FIM(j+1)
            for k=0: K_FIM(j+1)
                numAtom=numAtom+1;
            end
        end
    end
    dicParm = zeros(numAtom,3);
    indAtom=0;
    for j=0:L
        s = a^j;
        for p=0: P_FIM(j+1)
            u = p*(a^(j+1-l));
            for k=0: K_FIM(j+1)
                xi = k*(2*pi/N);
                indAtom=indAtom+1;
                dicParm(indAtom,:) = [s u xi];
            end
        end
    end
elseif strcmp(options,'allTimeShift')
    for j=0:L
        P_FIM(j+1) = N-1;
        K_FIM(j+1) = ceil((a^(l-1+j))/2) - 1;
        for k=0: K_FIM(j+1)
            for p=0: P_FIM(j+1)            
                numAtom=numAtom+1;
            end
        end
    end
    dicParm = zeros(numAtom,3);
    indAtom=0;
    for j=0:L
        s = a^j;
        for k=0: K_FIM(j+1)
            xi = k*pi*a^(2-j-l);
            for p=0: P_FIM(j+1)
                u = p;
                indAtom=indAtom+1;
                dicParm(indAtom,:) = [s u xi];
            end
        end
    end
end


end