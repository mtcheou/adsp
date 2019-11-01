function [signal sbaheader] = loadSBAfile (fileName) 

fid  = fopen(fileName,'r');

% Sign. Type
string = fgets(fid);
sigtype = sscanf(string(14:end),'%d');
% Dict. Type
string = fgets(fid);
dictype = sscanf(string(14:end),'%d');
% No. Signals
string = fgets(fid);
numSignal = sscanf(string(14:end),'%d');
% Signal Size
string = fgets(fid);
signalSize = sscanf(string(14:end),'%d');
% Block Hop
string = fgets(fid);
blockHop = sscanf(string(14:end),'%d');
% Block Size
string = fgets(fid);
blockSize = sscanf(string(14:end),'%d');
% Samp. Freq
string = fgets(fid);
Fs = sscanf(string(14:end),'%d');

sbaheader.sigtype = sigtype;
sbaheader.dictype = dictype;
sbaheader.numSignal = numSignal;
sbaheader.signalSize = signalSize;
sbaheader.blockHop = blockHop;
sbaheader.blockSize = blockSize;
sbaheader.Fs = Fs;


iSignal = 0;
iBlock = 0;
nMaxSignal = 1;
%block = struct('norm',0,'structbook',[]);
while(1)
    string = fgets(fid);
    if (string==-1) 
        break;
    elseif (string(1:4)=='XXXX')
        iSignal = iSignal+1;
        string = fgets(fid);
        nSignal = sscanf(string(9:end),'%d');
        string = fgets(fid);
        normSignal = sscanf(string(9:end),'%f');
        signal(nSignal).norm = normSignal;
    elseif (string(1:5)=='-----')
        iBlock = iBlock+1;
        string = fgets(fid);
        nBlock = sscanf(string(7:end),'%d');
        string = fgets(fid);
        normBlock = sscanf(string(7:end),'%f');
        string = fgets(fid);
        signal(nSignal).block(nBlock).norm = normBlock;
    elseif (string(1:5)=='#####')
        disp(['Block ' num2str(nBlock) ' with null samples' ])
    else
        sb = sscanf(string,'%f');
        if ( (sb(1)~=99999) && (sb(1)~=88888) && (sb(1)~=77777))
            signal(nSignal).block(nBlock).structbook(sb(1)).coef = sb(2);
            signal(nSignal).block(nBlock).structbook(sb(1)).rho = sb(3);
            signal(nSignal).block(nBlock).structbook(sb(1)).xi = sb(4);
            signal(nSignal).block(nBlock).structbook(sb(1)).phase = sb(5);
            signal(nSignal).block(nBlock).structbook(sb(1)).ti = sb(6);
            signal(nSignal).block(nBlock).structbook(sb(1)).tf = sb(7);
            signal(nSignal).block(nBlock).structbook(sb(1)).prevatom = sb(8);
        end
    end
end
   
    
fclose(fid);

