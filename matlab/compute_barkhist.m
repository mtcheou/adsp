function [cell_FLTLcount cell_FLTLenergy band_bark_all energy_all TLincrease_continuity FLincrease_continuity] = compute_barkhist(signal_band, band_bark_thresh,N)

band_bark_all = [];
energy_all=[];

for k1 = 1:length(signal_band)
    vec_FLTLcount=zeros(length(signal_band(k1).block),2);
    vec_FLTLenergy=zeros(length(signal_band(k1).block),2);
    TLincrease_continuity = 0;
    FLincrease_continuity = 0;
    for k2 = 1: length(signal_band(k1).block)
        FLcount=0;
        TLcount=0;
        FLenergy=0;
        TLenergy=0;
        for k3 = 1:length(signal_band(k1).block(k2).structbook)
            band_bark_all = [band_bark_all signal_band(k1).block(k2).structbook(k3).band_bark];
            energy_all = [energy_all (signal_band(k1).block(k2).structbook(k3).coef)^2];
            if (signal_band(k1).block(k2).structbook(k3).band_bark>band_bark_thresh)
                TLcount = TLcount+1;
                TLenergy = TLenergy + (signal_band(k1).block(k2).structbook(k3).coef)^2;
                if ( (signal_band(k1).block(k2).structbook(k3).rho<0) & ...
                     (signal_band(k1).block(k2).structbook(k3).tf==N) )
                     TLincrease_continuity = TLincrease_continuity + 1;
                end
            else
                FLcount = FLcount+1;
                FLenergy = FLenergy + (signal_band(k1).block(k2).structbook(k3).coef)^2;
                if ( (signal_band(k1).block(k2).structbook(k3).rho<0) & ...
                     (signal_band(k1).block(k2).structbook(k3).tf==N) )
                     FLincrease_continuity = FLincrease_continuity + 1;
                end
            end
        end
        vec_FLTLcount(k2,:) = [FLcount TLcount];
        vec_FLTLenergy(k2,:) = [FLenergy TLenergy];
    end
    cell_FLTLcount{k1} = vec_FLTLcount;
    cell_FLTLenergy{k1} = vec_FLTLenergy;
end