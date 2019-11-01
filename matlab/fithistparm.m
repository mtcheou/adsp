% fithistparm

fid = fopen('files.txt');
fileName = textscan(fid, '%s\n');
fileName =fileName{:};
fclose(fid);

[r c] = size(fileName);

parm_amp = [];
parm_rho = [];

for k = 1:r
    
    fname = strcat( '../', char(fileName(k)), 'strbookgrp.sbb');
    fname_header = strcat( '../', char(fileName(k)), 'strbookgrp_header.sbb');
    disp(['File: ' fname '; fname_header: ' fname_header]);
    
    [sbbHeader blockNorm structBook] = loadFileSBBGrp(fname_header, fname);
    
    for nsig = 1 : sbbHeader.numSignal
        for nblock = 1 : sbbHeader.numBlock
            sb = structBook{nsig,nblock};
            parm_amp = [parm_amp; sb(:,2)];
            parm_rho = [parm_rho; sb(:,3)];
        end
    end
    
end


figure, hist(parm_amp,1000);

figure, hist(parm_rho,1000);


