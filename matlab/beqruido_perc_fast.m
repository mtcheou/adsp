function [wi wf] = beqruido_perc_fast(x,perc,Fs)
% noise equivalent bandwidth (one-side)

N=length(x);
w = 0:(2*pi)/N: (N-1)*((2*pi)/N);
XX = abs(fft(x));
powtot = sum(XX(1:(N/2)+1).^2);
%disp(['powtot:' num2str(powtot)])
% calculate wi
k=1;
acum = 0;
while(1)
    acum = acum + sum(XX(k).^2);
%    disp(['Acum:' num2str(acum) '  ' num2str(perc*(powtot/2))])
    if (acum>perc*(powtot/2)) 
        break;
    end
    k=k+1;
end
%disp(['Acum left:' num2str(acum)])
ki=k;
wi = w(k);
% calculate wf
k=(N/2)+1;
acum=0;
while(1)
    acum = acum + sum(XX(k)^2);
%    disp(['Acum:' num2str(acum) '  ' num2str(perc*(powtot/2))])
    if (acum>perc*(powtot/2)) 
        break;
    end
    k=k-1;
end
%disp(['Acum right:' num2str(acum)])
kf=k;
wf = w(k);
