function x = synth_sb(signal,sbaheader,sbind)

blockSize = sbaheader.blockSize
blockHop = sbaheader.blockHop
signalSize = sbaheader.signalSize



x = zeros(length(signal),ceil(signalSize/blockHop)*blockHop); 

for k1 = 1:length(signal)
    for k2 = 1: length(signal(k1).block)
        sb = get_sbmatrix(k1,k2,signal,sbind);
        [r c] = size(sb);
        for k3 = 1: r
            realAtom = genexp(sb(k3,2:6),blockSize);
            
            %disp([ num2str(((k2-1)*blockHop)+1) ' '  num2str((k2-1)*blockHop+(blockSize))])
            x(k1,((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) =             ...
                       x(k1,((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) +  ...
                       (signal(k1).norm  * sb(k3,1) * realAtom); %* signal(k1).block(k2).norm
        end
    end
end