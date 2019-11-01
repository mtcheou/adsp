function [sbbHeader blockNorm structBook] = loadFileSBBGrp(fileHeader, fileSBB)

%
% INPUT -   fileHeader
%           fileSBB
%
% OUTPUT -  sbbHeader 
%           blockNorm 
%           structBook
%


%fileHeader = 'header.sbb';
%fileSBB = 'pianoA3_b1-11.sbb';

fid = fopen(fileHeader, 'r');

% Read header
sbbHeader.type = fread(fid, 1, 'int');

sbbHeader.dicType = fread(fid, 1, 'int');

sbbHeader.numSignal = fread(fid, 1, 'int');

sbbHeader.signalSize = fread(fid, 1, 'int');

sbbHeader.blockHop = fread(fid, 1, 'int');

sbbHeader.blockSize = fread(fid, 1, 'int');

sbbHeader.Fs =fread(fid, 1, 'double');

sbbHeader.numBlock = ceil(sbbHeader.signalSize/sbbHeader.blockHop);

fclose(fid);


% Read SBB file

fid = fopen(fileSBB, 'r');

disp('- Loading SBB File');



blockNorm = zeros(sbbHeader.numSignal,sbbHeader.numBlock);

sbbHeader.norm = zeros(1,sbbHeader.numSignal);


for nSignal=1:sbbHeader.numSignal
    disp(['nSignal: ' num2str(nSignal) ]);
    sbbHeader.norm(nSignal) = fread(fid, 1, 'double');    
    disp(['normSignal: ' num2str(sbbHeader.norm(nSignal))]);
    for nBlock =1: sbbHeader.numBlock
        disp(['nBlock: ' num2str(nBlock) ]);

        blockNorm(nSignal,nBlock) = fread(fid, 1, 'double');
        disp(['blockNorm: ' num2str(blockNorm(nSignal,nBlock))])
        
        nSBElement = fread(fid, 1, 'int');
        disp(['blockSignal: ' num2str(blockNorm(nSignal,nBlock))])
        
        sb = zeros(nSBElement,7);
        
        for k =1: nSBElement
            sb(k,1) = nSBElement;
            sb(k,2) = fread(fid, 1, 'double'); %   innerProduct
            sb(k,3) = fread(fid, 1, 'double'); %   rho
            sb(k,4) = fread(fid, 1, 'double'); %   xi
            sb(k,5) = fread(fid, 1, 'double'); %   phase
            sb(k,6) = fread(fid, 1, 'int'); % a 
            sb(k,7) = fread(fid, 1, 'int'); % b
            dummy = fread(fid, 1, 'int'); % nextAtom;
            dummy = fread(fid, 1, 'int'); % prevAtom;
            dummy = fread(fid, 1, 'int'); % origAtomIndex;
            dummy = fread(fid, 1, 'int'); % lack - folga;  
        end
        structBook{nSignal,nBlock}= sb;
    end
end

fclose(fid);

disp('- Loading has finished!!!');

