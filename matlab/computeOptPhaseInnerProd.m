function [opt_phi innerProd] = computeOptPhaseInnerProd(innerProd_xP,innerProd_xQ,innerProd_PQ,normQ2,normP2,xi)

% Input:
%       x - sinal de entrada
%       P - parte real do átomo complexo
%       Q - parte imaginária do átomo complexo
%       xi - freq. em radianos
%


a1 = innerProd_xP*normQ2 - innerProd_xQ*innerProd_PQ;
b1 = innerProd_xQ*normP2 - innerProd_xP*innerProd_PQ;

if (xi==0)
    opt_phi = 0;
    innerProd = innerProd_xP/sqrt(normP2);
elseif (xi~=0) && (a1==0)
    opt_phi = pi/2;
    innerProd = -innerProd_xQ/sqrt(norm(Q));
elseif (xi~=0) && (a1~=0)
    opt_phi = atan(-b1/a1);
    innerProd = (innerProd_xP*a1+innerProd_xQ*b1)/sqrt((a1^2)*normP2 + 2*a1*b1*innerProd_PQ + (b1^2)*normQ2);
end

end