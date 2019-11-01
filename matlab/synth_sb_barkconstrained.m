function [xFL xTL] = synth_sb_barkconstrained(signal,sbaheader,band_bark_thresh)

blockSize = sbaheader.blockSize
blockHop = sbaheader.blockHop
signalSize = sbaheader.signalSize

xFL_aux = zeros(length(signal),ceil(signalSize/blockHop)*blockHop); 
xTL_aux = zeros(length(signal),ceil(signalSize/blockHop)*blockHop);

for k1 = 1:length(signal)
    for k2 = 1: length(signal(k1).block)
        for k3 = 1:length(signal(k1).block(k2).structbook)
            p(1) = signal(k1).block(k2).structbook(k3).rho;
            p(2) = signal(k1).block(k2).structbook(k3).xi;
            p(3) = signal(k1).block(k2).structbook(k3).phase;
            p(4) = signal(k1).block(k2).structbook(k3).ti;
            p(5) = signal(k1).block(k2).structbook(k3).tf;                   
            coef = signal(k1).block(k2).structbook(k3).coef;
			realAtom = genexp(p,blockSize);    
            if (signal(k1).block(k2).structbook(k3).band_bark<band_bark_thresh)                        
            	xFL_aux(k1,((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) =          ...
                       	xFL_aux(k1,((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) +  ...
                       	(signal(k1).norm  * coef * realAtom);
            else
                xTL_aux(k1,((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) =          ...
                       	xTL_aux(k1,((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) +  ...
                       	(signal(k1).norm  * coef * realAtom);
            end
        end
    end
end

xFL = xFL_aux(1:signalSize);
xTL = xTL_aux(1:signalSize);