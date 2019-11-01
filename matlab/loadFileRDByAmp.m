function [header parm_range q_rho q_phi statFreqTable statFreqTableQ DfTable] = loadFileRDByAmp (fileName)


%
%   FUNCTION loadFileRDByAmp
%
%   INPUT - fileName
%
%
%   OUTPUT - header = [ numSignal 
%                       numBlock 
%                       signalSize 
%                       blockHop 
%                       blockSize
%                       maxNbitAmp 
%                       minNbitAmp 
%                       N_qrho 
%                       N_qphi];
%
%            parm_range(numSignal,numBlock) = struct('min_amp',[],'max_amp',[], ...
%                                                    'min_rho',[],'max_rho',[], ...
%                                                    'min_phase',[],'max_phase',[]);
%            q_rho
%
%            q_phi
%
%            statFreqTable= zeros((maxNbitAmp-minNbitAmp+1),maxNbitAmp);      
%
%            statFreqTableQ= zeros((maxNbitAmp-minNbitAmp+1),maxNbitAmp);
%
%            DfTable = zeros((maxNbitAmp-minNbitAmp+1),maxNbitAmp,N_qrho,N_qphi);
%


fid = fopen(fileName, 'r');

%fid = fopen('./DFTABLE_byamp/dftable_byamp_nbamp5-5.bin', 'r');


% Read header
header = struct('numSignal',[],...
                'numBlock',[], ...
                'signalSize',[],...
                'blockHop',[], ...
                'blockSize',[], ...
                'maxNbitAmp',[], ...
                'minNbitAmp',[], ...
                'Nqrho',[], ...
                'Nqphi',[]);

numSignal = fread(fid, 1, 'int');
numBlock = fread(fid, 1, 'int');
signalSize = fread(fid, 1, 'int');
blockHop = fread(fid, 1, 'int');
blockSize =fread(fid, 1, 'int');
maxNbitAmp =fread(fid, 1, 'int');
minNbitAmp =fread(fid, 1, 'int');
N_qrho =fread(fid, 1, 'int');
N_qphi = fread(fid, 1, 'int');

header.numSignal = numSignal;
header.numBlock = numBlock;
header.signalSize = signalSize;
header.blockHop = blockHop;
header.blockSize = blockSize;
header.maxNbitAmp = maxNbitAmp;
header.minNbitAmp = minNbitAmp;
header.Nqrho = N_qrho;
header.Nqphi = N_qphi;

% Read parameters ranges

parm_range(numSignal,numBlock) = struct('min_amp',[],'max_amp',[], ...
                                        'min_rho',[],'max_rho',[], ...
                                        'min_phase',[],'max_phase',[]);

for k1= 1:numSignal
    for k2= 1:numBlock
        parm_range(k1,k2).min_amp = fread(fid, 1, 'double');
        parm_range(k1,k2).max_amp = fread(fid, 1, 'double');
        parm_range(k1,k2).min_rho = fread(fid, 1, 'double');
        parm_range(k1,k2).max_rho = fread(fid, 1, 'double');
        parm_range(k1,k2).min_phase = fread(fid, 1, 'double');
        parm_range(k1,k2).max_phase = fread(fid, 1, 'double');
    end
end


% Read rho and phase quatization steps
q_rho = zeros(1,N_qrho);
for k1 = 1:N_qrho
    q_rho(k1)= fread(fid, 1, 'double');
end

q_phi = zeros(1,N_qphi);
for k1 = 1:N_qphi
    q_phi(k1)= fread(fid, 1, 'double');
end


% Read frequency table of non quantized atoms

statFreqTable= zeros((maxNbitAmp-minNbitAmp+1),maxNbitAmp);

for k1 =minNbitAmp:maxNbitAmp
    for k2=1:k1
        statFreqTable(k1-minNbitAmp+1,k2) = fread(fid, 1, 'int');   
    end
end


statFreqTableQ= zeros((maxNbitAmp-minNbitAmp+1),maxNbitAmp);

for k1 =minNbitAmp:maxNbitAmp
    for k2=1:k1
        statFreqTableQ(k1-minNbitAmp+1,k2) = fread(fid, 1, 'int');
    end
end

% Read the mean atom distortion atom for each amplitude range

% DfTable = zeros(numSignal,numBlock,(maxNbitAmp-minNbitAmp+1),maxNbitAmp,N_qrho,N_qphi);
% for k1 = 1:numSignal
%     for k2=1: numBlock
%         for k3 =minNbitAmp:maxNbitAmp
%             for k4=1:k3
%                 for k5 =1:N_qrho
%                     for k6=1:N_qphi
%                         DfTable(k1,k2,k3-minNbitAmp+1,k4,k5,k6) = fread(fid, 1, 'double');
%                     end
%                 end
%             end
%         end
%     end
% end

DfTable = zeros((maxNbitAmp-minNbitAmp+1),maxNbitAmp,N_qrho,N_qphi);
for k3 =minNbitAmp:maxNbitAmp
    for k4=1:k3
        for k5 =1:N_qrho
            for k6=1:N_qphi
                DfTable(k3-minNbitAmp+1,k4,k5,k6) = fread(fid, 1, 'double');
            end
        end
    end
end
 

fclose(fid);