% decompGabor
N = 32;

dicGamma = getDicParmDiscrete('std',N);

[row ,col] = size(dicGamma); 
numGamma = row;

dicGaborComplex = zeros(numGamma,N);

for k = 1:numGamma
    dicGaborComplex(k,:) = genGaborComplex( dicGamma(k,:),N);
end

x = genGaborReal([32 0 0.589],pi/8,N); 
%x = randn(1,N);
residue=x;

innerProd = zeros(1,numGamma);
opt_phi = zeros(1,numGamma);

NIter = 5;
chosenParm = zeros(NIter,5);
for kIter = 1: NIter
    for k = 1:numGamma
        Pgamma = real(dicGaborComplex(k,:));
        Qgamma = imag(dicGaborComplex(k,:));
        innerProd_xP = residue*Pgamma';
        innerProd_xQ = residue*Qgamma';
        innerProd_PQ = Pgamma*Qgamma';
        normP2 = norm(Pgamma)^2;
        normQ2 = Qgamma*Qgamma';
        xi = dicGamma(k,3);
        [opt_phi(k),innerProd(k)] = computeOptPhaseInnerProd(innerProd_xP,innerProd_xQ,innerProd_PQ,normQ2,normP2,xi);
    end
    innerProdMax = max(innerProd);
    indAtomMax = find(innerProd==innerProdMax);

    chosenAlpha = innerProdMax;
    chosenGamma = dicGamma(indAtomMax,:);
    chosenOptPhi = opt_phi(indAtomMax);
    chosenAtom = genGaborReal(chosenGamma,chosenOptPhi,N);
    chosenParm(kIter,:) = [chosenAlpha chosenGamma chosenOptPhi]; 
    residue = residue - chosenAlpha*chosenAtom;
end

