function evaldecomp_contsb_coeforder(signal,sbaheader,contAtom,origsignal,nContAtom)

blockSize = sbaheader.blockSize;
blockHop = sbaheader.blockHop;
signalSize = sbaheader.signalSize;

for k =1:length(contAtom)
    sb = contAtom{k};
    norm2(k) = (sum(sb(:,3)))^2;
end

[norm2_y ind] = sort(norm2,'descend');

count=0;
count2=1;

%Define the area to be recorded 
%hf=figure; 
%rect = get(hf,'Position'); 
%rect(1:2) = [0 0]; 
figure;
res = origsignal;
while (count2<=nContAtom)
    if (count2>length(contAtom)) 
        break;
    end
    x = zeros(1,ceil(signalSize/blockHop)*blockHop);
    sb = contAtom{ind(count2)};
    [r c] = size(sb);
    for k = 1:r
        realAtom = genexp(sb(k,4:8),blockSize);
        k2 = sb(k,1);
        x(((k2-1)*blockHop)+1:(k2-1)*blockHop+(blockSize)) = ...
                   (signal.norm  * sb(k,3) * realAtom);
    end
    plot(res);
    hold on;
    plot(x,'r--');
    hold off;
    res = res - x(1:signalSize)';
    title([' Continued Atom: ' num2str(count2) ' -- Norm_residue: ' num2str(norm(res))])
    % Generate and record the frames 
    %M(:,count2) = getframe(hf,rect); 
    pause;
    count = count + r;
    count2 = count2 + 1;
end

%save evaldecomp.mat M

%movie2avi(M, 'evaldecomp.avi', 'compression', 'None');




