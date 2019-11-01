function g_gamma = genGaborComplex(parm,N)


s = parm(1);
u = parm(2);
xi = parm(3);

g_gamma = zeros(1,N);

if (s==1)
    g_gamma(u+1) = 1 * exp(i*xi*(u));
elseif (s==N)
    for n= 0:N-1
        g_gamma(n+1) = exp(i*xi*n);
    end
else
    for n= 0:N-1
        g_gamma(n+1) = (2^(1/4)) * exp(-pi*((n-u)^2)) * exp(i*xi*n);
    end
end

g_gamma = g_gamma / norm(g_gamma);