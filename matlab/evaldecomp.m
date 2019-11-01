function evaldecomp (fileName, origSignal, signal, block) 

residue = origSignal;
z = origSignal/norm(origSignal);

fid  = fopen(fileName,'r');

% Sign. Type
string = fgets(fid);
sigtype = sscanf(string(14:end),'%d')
% Dict. Type
string = fgets(fid);
dictype = sscanf(string(14:end),'%d')
% No. Signals
string = fgets(fid);
nSignal = sscanf(string(14:end),'%d')
% Signal Size
string = fgets(fid);
signalSize = sscanf(string(14:end),'%d')
% Block Hop
string = fgets(fid);
blockHop = sscanf(string(14:end),'%d')
% Block Size
string = fgets(fid);
blockSize = sscanf(string(14:end),'%d')
% Samp. Freq
string = fgets(fid);
Fs = sscanf(string(14:end),'%d')
% Init Block
string = fgets(fid);
initBlock = sscanf(string(14:end),'%d')
% Final Block
string = fgets(fid);
finalBlock = sscanf(string(14:end),'%d')

iSignal = 0;
iBlock = 0;
nMaxSignal = 1;
x = zeros(1,blockSize);
x = origSignal(((block-1)*blockHop)+1:(block-1)*blockHop+(blockSize));
figure;
while(1)
    string = fgets(fid)
    if (string==-1) 
        break;
    elseif (string(1:4)=='XXXX')
        iSignal = iSignal+1;
        string = fgets(fid);
        nSignal = sscanf(string(9:end),'%d');
        string = fgets(fid);
        normSignal = sscanf(string(9:end),'%f');
    elseif (string(1:5)=='-----')
        iBlock = iBlock+1;
        string = fgets(fid);
        nBlock = sscanf(string(7:end),'%d');
        string = fgets(fid);
        normBlock = sscanf(string(7:end),'%f');
        string = fgets(fid);
    elseif (string(1:5)=='#####')
        disp(['Block ' num2str(nBlock) ' with null samples' ])
    else
        sb = sscanf(string,'%f')
        if ( (sb(1)~=99999) && (sb(1)~=88888) && (sb(1)~=77777) && (block==nBlock) && (signal==nSignal))
            realAtom = genexp(sb(3:7),blockSize);
            plot(x);
            hold on;
            plot(normSignal*sb(2)*realAtom,'r--');
            hold off;
            x = x - normSignal*sb(2)*realAtom';
            title(['Block: ' num2str(nBlock) '-- Atom: ' num2str(sb(1)) ' -- Norm_res: ' num2str(norm(x))])
            %axes([0 10000 -1 1])
            pause;
        end
    end
end
   
iSignal
iBlock
    
fclose(fid);

