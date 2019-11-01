

fid = fopen('files_dftable_byamp_complete.txt');
fileName = textscan(fid, '%s');
fileName =fileName{:};
fclose(fid);

[r c] = size(fileName);

fname = strcat ( './', char(fileName(1)) );

[header parm_range q_rho q_phi statFreqTable statFreqTableQ DfTable] = loadFileRDByAmp (fname);

totalStatFreqTableQ = zeros(16,16);
for k2 = 1:16    
    for k3 = 1:k2
        meanDfTable{k2,k3} = zeros(header.Nqrho,header.Nqphi);
    end
end

for k = 1:r
    
    fname = strcat ( './',  char(fileName(k)) );
    disp(['File: ' fname]);
    [header parm_range q_rho q_phi statFreqTable statFreqTableQ DfTable] = loadFileRDByAmp (fname);
    disp(['minNbAmp: ' num2str(header.minNbitAmp) ]);
    disp(['maxNbAmp: ' num2str(header.maxNbitAmp) ]);
    for k2 = header.minNbitAmp:header.maxNbitAmp
        disp(['nbamp: ' num2str(k2) ]);
        for k3 = 1:k2
            disp(['iAmpRange: ' num2str(k3) ]);
            
            totalStatFreqTableQ(k2,k3) = totalStatFreqTableQ(k2,k3) + statFreqTableQ(k2-header.minNbitAmp+1,k3);
            disp(['statFreqTableQ: ' num2str(statFreqTableQ(k2-header.minNbitAmp+1,k3)) ]);
            disp(['totalStatFreqTableQ: ' num2str(totalStatFreqTableQ(k2,k3)) ]);
            
            if (statFreqTableQ(k2-header.minNbitAmp+1,k3)~=0)
                x = reshape(DfTable(k2-header.minNbitAmp+1,k3,:,:),header.Nqrho,header.Nqphi);
                meanDfTable{k2,k3} = meanDfTable{k2,k3} + statFreqTableQ(k2-header.minNbitAmp+1,k3)*x;
            end
        end
        
    end
    
end

for k2 = 1:16
    figure;
    delta = 16/( 2^(ceil(log2(k2))));
    %title(['NbAmp: ' num2str(k2)]);
    for k3 = 1:k2
        meanDfTable{k2,k3} = meanDfTable{k2,k3}/totalStatFreqTableQ(k2,k3) ;
        
        %meanDfTable_dB{k2,k3} = db(x,'power');
        
        plotplace = [1:delta] + (k3-1)*delta;
        x = meanDfTable{k2,k3};
        subplot(4,4,plotplace),imagesc(q_phi,q_rho,db(x,'power')); colorbar; 
        %colormap(flipud(gray));
        %colormap(gray);
        xlabel('qphi')
        ylabel([ 'Range:' num2str(k3) '- qrho'])
        %pause;
    end
end

fid = fopen('dftablebyamp.bin','wb');

fwrite(fid, header.Nqrho, 'int');
fwrite(fid, header.Nqphi, 'int');

N_qrho = header.Nqrho;
N_qphi = header.Nqphi;

for k =1:N_qrho
    fwrite(fid, q_rho(k), 'double');
end

for k=1:N_qphi
    fwrite(fid, q_phi(k), 'double');
end


for k2 =1:16
    for k3=1:k2
        x = meanDfTable{k2,k3};
        %imagesc(q_phi,q_rho,x); colorbar;
        for k5 =1:N_qrho
            for k6=1:N_qphi 
                fwrite(fid, x(k5,k6), 'double');
            end
        end
        %pause;
    end
end

fclose(fid);

% 
% meanDfTabInterp=zeros(N_qrho,N_qphi);
% A = zeros(N_qrho*N_qphi,6);
% for k2 =1:16
%     for k3=1:k2
%         meandftab = meanDfTable{k2,k3};
%         B = DCT2(meandftab);
%         
%         k=1;
%         for k5 =1:N_qrho
%             for k6=1:N_qphi 
%                 X = [  q_rho(k5)^2 ...
%                        q_rho(k5)*q_phi(k6) ...
%                        q_phi(k6)^2 ...
%                        q_rho(k5) ...
%                        q_phi(k6) ...
%                        1    ];
%                 A(k,:) = meandftab(k5,k6)*pinv(X);
%                 k=k+1;           
%             end
%         end
%         Amean = mean(A);
%         k=1;
%         for k5 =1:N_qrho
%             for k6=1:N_qphi 
%                 X = [  q_rho(k5)^2 ...
%                        q_rho(k5)*q_phi(k6) ...
%                        q_phi(k6)^2 ...
%                        q_rho(k5) ...
%                        q_phi(k6) ...
%                        1    ];
%                 %meanDfTabInterp(k5,k6) = X*Amean'; 
%                 meanDfTabInterp(k5,k6) = X*A(k,:)'; 
%                 k=k+1
%             end
%         end
%         meanDfTableInterp{k2,k3}=meanDfTabInterp;
%         subplot(2,1,1),imagesc(db(meanDfTable{k2,k3})); colorbar;
%         subplot(2,1,2),imagesc(db(meanDfTableInterp{k2,k3})); colorbar;
%         title(['NBAmp: ' num2str(k2) ' AmpRange: ' num2str(k3)])
%         pause;
%     end
% end
% 




            