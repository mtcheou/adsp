% compute_band_sb
function signal_band = compute_band_sb(signal,sbaheader)
% matrix columns : 

N=sbaheader.blockSize;
Nzeros=N*32;
Fs = sbaheader.Fs;
p = zeros(1,5);
x =zeros(1,N+Nzeros);
for k1 = 1:length(signal)
    for k2 = 1: length(signal(k1).block)
        for k3 = 1: length(signal(k1).block(k2).structbook)
            disp(['signal: ' num2str(k1) ' block: ' num2str(k2) ' SB element: ' num2str(k3)])
            p(1) = signal(k1).block(k2).structbook(k3).rho;
            p(2) = signal(k1).block(k2).structbook(k3).xi;
            p(3) = signal(k1).block(k2).structbook(k3).phase;
            p(4) = signal(k1).block(k2).structbook(k3).ti;
            p(5) = signal(k1).block(k2).structbook(k3).tf;
            % compute band
            x(1:N) = genexp(p,N);
            [wi wf] = beqruido_perc_fast(x,0.3,Fs);
            band_bark = (bark((Fs/(2*pi))*wf) - bark((Fs/(2*pi))*wi));
            % save in the same structure
            signal(k1).block(k2).structbook(k3).wi = wi;
            signal(k1).block(k2).structbook(k3).wf = wf;
            signal(k1).block(k2).structbook(k3).band_bark = band_bark;
        end
    end
end

signal_band = signal;