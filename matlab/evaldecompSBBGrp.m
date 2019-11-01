function evaldecompSBBGrp(fileHeader, fileSBB, origSignal, signal, block, type, coeforder)
%
%   evaldecompSBB(fileHeader, fileSBB, origSignal, signal, block, type, coeforder)
%
%       fileHeader
%       fileSBB
%       origSignal (length,signal)
%       signal
%       block
%       type - 'decomp': decomposition step by step
%            - 'resnorm': residue norm wrt to decomp. iteration
%            - 'snr': signal to noise ratio wrt to decomp. iteration
%            - 'coef' : coef wrt to decomp. iteration
%       coeforder - 1 considered; 0 not considered
%

[sbbHeader blockNorm structBook] = loadFileSBBGrp(fileHeader, fileSBB);


normSignal = norm(origSignal(:,signal));
z = origSignal(:,signal)/normSignal;


sb = structBook{signal,block};
[r c] = size(sb);

if coeforder
    sb = sortrows(sb,-2);
end

initBlock = ((block-1)*sbbHeader.blockHop)+1;
endBlock = (block-1)*sbbHeader.blockHop+(sbbHeader.blockSize);
if (endBlock > sbbHeader.signalSize)
    endBlock = sbbHeader.signalSize;
end

x = zeros(1,sbbHeader.blockSize);
residue = x;
x(1:endBlock-initBlock+1) = z(initBlock:endBlock);
residue(1:endBlock-initBlock+1) = origSignal(initBlock:endBlock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(type, 'decomp')
    for k = 1:r
        realAtom = genexp(sb(k,3:7),sbbHeader.blockSize);
        plot(residue);
        set(gca,'XTick',[1:4096:sbbHeader.signalSize],'XGrid','on')
        hold on;
        plot(sbbHeader.norm(signal)*sb(k,2)*realAtom,'r--');
        hold off;
        residue = residue - sbbHeader.norm(signal)*sb(k,2)*realAtom;
        title(['Signal: ' num2str(signal) '--Block: ' num2str(block) '-- Atom: ' num2str(sb(k,1)) ' -- Norm_res: ' num2str(norm(residue))])
        pause;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(type,'resnorm')
    resnorm = zeros(1,r+1);
    for k = 1:r
        resnorm(k) = norm(x);
        realAtom = genexp(sb(k,3:7),sbbHeader.blockSize);
        x = x - sb(k,2)*realAtom;
    end
    resnorm(r+1) = norm(x);
    figure, plot(0:r,resnorm);
    title(['Signal: ' num2str(signal) '-- Block: ' num2str(block)])
    xlabel('Iteration')
    ylabel('Residue norm');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(type, 'snr')
    snr_dB = zeros(1,r+1);
    for k = 1:r
        snr_dB(k) = 20*log10(blockNorm(signal,block)/norm(x));
        realAtom = genexp(sb(k,3:7),sbbHeader.blockSize);
        size(realAtom);
        x = x - sb(k,2)*realAtom;
    end
    snr_dB(r+1) = 20*log10(blockNorm(signal,block)/norm(x));
    figure, plot(0:r,snr_dB);
    title(['Signal: ' num2str(signal) '-- Block: ' num2str(block)])
    xlabel('Iteration')
    ylabel('SNR (dB)');
end

if strcmp(type, 'coef')
    figure, plot(1:r,sb(:,2));
    title(['Signal: ' num2str(signal) '-- Block: ' num2str(block)])
    xlabel('Iteration')
    ylabel('Amplitude');
end



