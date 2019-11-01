% plotLagrangian.m

function plotLagrangian(signal,sbaheader,nSignal,nBlock)


sb = zeros(length(signal(nSignal).block(nBlock).structbook),6);

for k =1: length(signal(nSignal).block(nBlock).structbook)
    sb(k,:) = [   signal(nSignal).block(nBlock).structbook(k).coef ...
                signal(nSignal).block(nBlock).structbook(k).rho ...
                signal(nSignal).block(nBlock).structbook(k).xi ...
                signal(nSignal).block(nBlock).structbook(k).phase ...
                signal(nSignal).block(nBlock).structbook(k).ti ...
                signal(nSignal).block(nBlock).structbook(k).tf];
end

[r c] = size(sb);

max_rho = max(abs(sb(:,2)))

max_phi = 2*pi

delta = 1/((2^3)-1);

qrho = 0:delta:1;
qphi = 0:delta:1;

mse = zeros(length(qrho),length(qphi));
y = zeros(length(qrho),length(qphi));


lambda=0.1;

for k1 = 1:length(qrho)
    for k2 = 1:length(qphi)
        disp(['qrho: ' num2str(k1) 'qphi: ' num2str(k2)])
        for k= 1: r
            sbq = sb(k,:);

            if (sb(k,2)>=0)
                sbq(2) = linearquant(abs(sb(k,2)),max_rho*qrho(k1));
            else
                sbq(2) = - linearquant(abs(sb(k,2)),max_rho*qrho(k1));
            end

            sbq(4) = linearquant(sb(k,4),max_phi*qphi(k2));
            
            [origAtom normr] = genexp(sb(k,2:end),sbaheader.blockSize);
            [quantAtom normr] = genexp(sbq(2:end),sbaheader.blockSize);
            
%             figure, plot(origAtom)
%             hold on
%             plot(quantAtom,'r')
%             pause;
%             close;
            
            mse(k1,k2) = mse(k1,k2) + ((norm(origAtom-quantAtom))^2)/sbaheader.blockSize;
        end
        mse(k1,k2) = mse(k1,k2)/r;
        y(k1,k2) = mse(k1,k2)/r  + lambda *(log2(1/qrho(k1))+ log2(1/qphi(k2))); 
    end
end

save mse.mat mse qrho qphi y lambda

figure,contour(qrho,qphi,mse); 

figure,contour(qrho,qphi,mse); 
    
   
   






