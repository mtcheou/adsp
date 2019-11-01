function g_gamma = genGaborReal(parm,phi,N)

s = parm(1);
u = parm(2);
xi = parm(3);

g_gamma = zeros(1,N);

if (s==1)
    g_gamma(u+1) = 1 ;
elseif (s==N)
    for n= 0:N-1
        g_gamma(n+1) = cos(xi*n + phi);
    end
else
    for n= 0:N-1
        g_gamma(n+1) = (2^(1/4)) * exp(-pi*((n-u)^2)) * cos(xi*n + phi);
    end
end

g_gamma = g_gamma / norm(g_gamma);