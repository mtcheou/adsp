function x = plot_blockatom (fileName, origSignal, nAtom) 

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

iSignal = 0;
iBlock = 0;
nMaxSignal = 1;
x = zeros(nAtom,signalSize); 
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
        sb = sscanf(string,'%f');
        if ( (sb(1)~=9999) && (sb(1)~=8888) && (sb(1)~=7777) && (sb(1)<=nAtom))
%             normSignal
%             normBlock
%             sb(2)'
%             sb(3:7)'
            realAtom = genexp(sb(3:7),blockSize);
            x(sb(1),((nBlock-1)*blockHop)+1:(nBlock-1)*blockHop+(blockSize)) =      ...
                normSignal*sb(2)*realAtom;
%            pause;
        end
    end
end
   
iSignal
iBlock
    
fclose(fid);

