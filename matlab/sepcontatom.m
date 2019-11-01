function contAtom = sepcontatom (signal,sbaheader)

blockSize = sbaheader.blockSize;
blockHop = sbaheader.blockHop;
signalSize = sbaheader.signalSize;

numBlock = ceil(signalSize/blockSize);
flagcont = zeros(2000,numBlock);
nContAtom=1;
for k1 = 1:length(signal)
    for k2 = length(signal(k1).block):-1:1
        for k3 = 1: length(signal(k1).block(k2).structbook)
            if (flagcont(k3,k2)==0)
                flagcont(k3,k2)=1;
                contsb = [ ]; 
                k2aux = k2;                
                prevatom = k3;
                while(k2aux>0)
                    contsb = [ contsb ;  ...
                                 k2aux ...
                                 signal(k1).block(k2aux).norm ...
                                 signal(k1).block(k2aux).structbook(prevatom).coef ...
                                 signal(k1).block(k2aux).structbook(prevatom).rho ...
                                 signal(k1).block(k2aux).structbook(prevatom).xi ...
                                 signal(k1).block(k2aux).structbook(prevatom).phase ...
                                 signal(k1).block(k2aux).structbook(prevatom).ti ...
                                 signal(k1).block(k2aux).structbook(prevatom).tf];
                    prevatom = signal(k1).block(k2aux).structbook(prevatom).prevatom; 
                    if (prevatom==0) 
                        break;
                    end
                    k2aux = k2aux -1;
                    flagcont(prevatom,k2aux)=1; 
                end
                contAtom{nContAtom} = contsb;
                nContAtom=nContAtom+1;
            end
        end
    end
end


