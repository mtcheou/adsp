function x = synth_sb_coeforder(signal,sbaheader,nAtom)

blockSize = sbaheader.blockSize
blockHop = sbaheader.blockHop
signalSize = sbaheader.signalSize



x = zeros(length(signal),ceil(signalSize/blockHop)*blockHop); 

for k1 = 1:length(signal)
    for k2 = 1: length(signal(k1).block)
        sb = zeros(length(signal(k1).block(k2).structbook),6);
        for k3 = 1: length(signal(k1).block(k2).structbook)
            sb(k3,:) = [ signal(k1).block(k2).structbook(k3).coef ...
                 signal(k1).block(k2).structbook(k3).rho ...
                 signal(k1).block(k2).structbook(k3).xi ...
                 signal(k1).block(k2).structbook(k3).phase ...
                 signal(k1).block(k2).structbook(k3).ti ...
                 signal(k1).block(k2).structbook(k3).tf];
        end
        %figure,plot(sb(:,1));
        sb = flipud(sortrows(sb,1));
        %hold on;
        %plot(sb(:,1),'r');
        if (~isempty(sb))
            for k3 = 1: nAtom
                realAtom = genexp(sb(k3,2:6),blockSize);
                
                %disp([ num2str(((k2-1)*blockHop)+1) ' '  num2str((k2-1)*blockHop+(blockSize))])
                x(k1,((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) =             ...
                           x(k1,((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) +  ...
                           (signal(k1).norm  * sb(k3,1) * realAtom); %* signal(k1).block(k2).norm
            end
         end
    end
end