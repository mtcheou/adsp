function [sbbHeader blockNorm structBook] = loadFileSBB(fileHeader, fileSBB)

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

sb = [];

blockNorm = zeros(sbbHeader.numSignal,sbbHeader.numBlock);

sbbHeader.norm = zeros(1,sbbHeader.numSignal);

while(1)
    nSignal = fread(fid, 1, 'int');
    if (nSignal==77777) 
        disp('- Loading complete!!!');
        break;
    end
    sbbHeader.norm(nSignal) = fread(fid, 1, 'double');
    while(1)
        nBlock = fread(fid, 1, 'int');
        if (nBlock==88888) 
            disp('- Signal loading complete!!!');
            break;
        end
        blockNorm(nSignal,nBlock) = fread(fid, 1, 'double');
        while(1)
            nSBElement = fread(fid, 1, 'int');
            if (nSBElement==99999)
                structBook{nSignal,nBlock}= sb;
                sb=[];
                disp('- Block loading complete!!!');
                break;
            end
            %disp(['nSBElement: ' num2str(nSBElement)]);
            sb(nSBElement,1) = nSBElement;
            sb(nSBElement,2) = fread(fid, 1, 'double'); %   innerProduct
            sb(nSBElement,3) = fread(fid, 1, 'double'); %   rho
            sb(nSBElement,4) = fread(fid, 1, 'double'); %   xi
            sb(nSBElement,5) = fread(fid, 1, 'double'); %   phase
            sb(nSBElement,6) = fread(fid, 1, 'int'); % a 
            sb(nSBElement,7) = fread(fid, 1, 'int'); % b
            dummy = fread(fid, 1, 'int'); % nextAtom;
            dummy = fread(fid, 1, 'int'); % prevAtom;
            dummy = fread(fid, 1, 'int'); % origAtomIndex;
            dummy = fread(fid, 1, 'int'); % lack - folga;  
        end
    end
end

fclose(fid);

disp('- Loading has finished!!!');

