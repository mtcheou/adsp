function x = synth_contsb_coeforder(signal,sbaheader,contAtom,nAtom)

blockSize = sbaheader.blockSize;
blockHop = sbaheader.blockHop;
signalSize = sbaheader.signalSize;

for k =1:length(contAtom)
    sb = contAtom{k};
    norm2(k) = (sum(sb(:,3)))^2;
end

[norm2_y ind] = sort(norm2,'descend');

x = zeros(1,ceil(signalSize/blockHop)*blockHop); 

count=0;
count2=1;
while (count<nAtom-1)
    if (count2>length(contAtom)) 
        break;
    end
    sb = contAtom{ind(count2)};
    [r c] = size(sb);
    for k = 1:r
        realAtom = genexp(sb(k,4:8),blockSize);
        k2 = sb(k,1);
        x(((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) =             ...
                   x(((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) +  ...
                   (signal.norm  * sb(k,3) * realAtom);
    end
    count = count + r;
    count2 = count2 + 1;
end
count
count2