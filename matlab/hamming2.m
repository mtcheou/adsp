function w = hamming2(L,S)
a = .54;
b = -.46;

w = zeros(L,1);
K = (2*sqrt(S/L))/sqrt(4*a^2+2*b^2);

for n = 0:L-1
    w(n+1) = K * (a + b*cos(2*pi*n/L+pi/L));
end
    
end