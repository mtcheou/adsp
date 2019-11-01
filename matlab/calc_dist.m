%calc_dist


dphi =2^(-5);
drho = 2^(-5);

load sblinked_load.out

sb = sblinked_load;

%sb = [0.000407520000000000 0.000435920000000000 0.0107665400000000 5.38723239000000 28241 36863];

[r c] =size(sb);

Dfcomp = zeros(1,r);

qmat = zeros(r,4);

for k= 1:r
    rho = sb(k,2);
    if (rho<0) 
        rho = -rho;
    end
    w = sb(k,3);
    phi = sb(k,4);
    dt = sb(k,6) - sb(k,5);
    signalSize = ceil(dt/4096)*4096;
    [realAtom normr(k)] = genexp([rho w phi 0 dt],dt+1);
    
    rhoq = linearquant(rho,drho);
    phiq = linearquant(phi,dphi);
    
    qmat(k,:) = [rho rhoq phi phiq];
    
    [realAtomq normrq(k)] = genexp([rhoq w phiq 0 dt],dt+1);
    
    Dfcomp(k) = sum((normr(k)*realAtom-normrq(k)*realAtomq).^2 )/(normr(k)*normr(k));
end
    

figure ,plot(Dfcomp);