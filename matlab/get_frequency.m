function [ count_block sb_freqlike sb_timelike] = get_frequency(signal,sbaheader)
%
% get_frequency: separates frequency-like atoms from time-like atoms 
%
% [count_block sb_freqlike sb_timelike] = get_frequency(signal,sbaheader)
%
%

load freq_like_table.mat

Fs = sbaheader.Fs;
N = sbaheader.blockSize;

for k1 = 1:length(signal)
    count_block = zeros(length(signal(k1).block),4);
    for k2 = 1: length(signal(k1).block)
        count_freqlike = 0;
        count_timelike = 0;
        energy_freqlike =0;
        energy_timelike =0;
        sbindex_freqlike = [];
        sbindex_timelike = [];
        for k3 = 1: length(signal(k1).block(k2).structbook)
            if (signal(k1).block(k2).structbook(k3).rho ==0)
                s =N;
            else   
                s = 2^(min ([floor(log2(abs(1/signal(k1).block(k2).structbook(k3).rho))) log2(N)]));
            end
            mat = mat_result(find(mat_result(:,1)==s),:);
            [r c ] = size(mat);
            freq = (Fs/(2*pi))*signal(k1).block(k2).structbook(k3).xi;
            index = find(abs(mat(:,2) - freq*ones(r,1))==min(abs( mat(:,2) - freq*ones(r,1))));
            delta_tau = signal(k1).block(k2).structbook(k3).tf - ...
                        signal(k1).block(k2).structbook(k3).ti;
            if (isempty(mat))
                count_timelike = count_timelike + 1;
                energy_timelike = energy_timelike + signal(k1).block(k2).structbook(k3).coef;
                sbindex_timelike = [sbindex_timelike; k3];
            else
                if ( (freq<mat(1,2)) )
                    count_timelike = count_timelike + 1;
                    energy_timelike = energy_timelike + signal(k1).block(k2).structbook(k3).coef;
                    sbindex_timelike = [sbindex_timelike; k3];
                else 
                    if (delta_tau > (N - mat(index(1),3)))
                        count_freqlike = count_freqlike + 1;
                        energy_freqlike = energy_freqlike + signal(k1).block(k2).structbook(k3).coef;
                        sbindex_freqlike = [sbindex_freqlike; k3];
                    else
                        count_timelike = count_timelike + 1;
                        energy_timelike = energy_timelike + signal(k1).block(k2).structbook(k3).coef;
                        sbindex_timelike = [sbindex_timelike; k3];
                    end
                end
            end
        end
        count_block(k2,:) = [count_freqlike energy_freqlike count_timelike energy_timelike];
        sb_freqlike{k1,k2} = sbindex_freqlike; 
        sb_timelike{k1,k2} = sbindex_timelike;
    end
end