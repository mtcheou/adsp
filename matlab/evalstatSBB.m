function evalstatSBB(fileHeader, fileSBB, signal, block, type, parm, nbin, dist)
%
%   evalstatSBB(fileHeader, fileSBB, signal, block, type, parm, nbin)
%       fileHeader
%       fileSBB
%       signal
%       block
%       type - 'hist'
%            - 'blocknorm'
%            - 'parmrange'
%       parm - 'coef'
%            - 'decay'  
%            - 'freq'  
%            - 'phase'  
%            - 'timesupport'  
%       nbin- number of histogram bins
%       dist - ditribution function
%                 'beta'
%                 'birnbaumsaunders'
%                 'exponential'
%                 'extreme value' or ev'
%                 'gamma'
%                 'generalized extreme value' or 'gev'
%                 'generalized pareto' or 'gp'
%                 'inversegaussian'
%                 'logistic'
%                 'loglogistic'
%                 'lognormal'
%                 'nakagami'
%                 'negative binomial' or 'nbin'
%                 'normal' (default)
%                 'poisson'
%                 'rayleigh'
%                 'rician'
%                 'tlocationscale'
%                 'weibull' or 'wbl'
%


[sbbHeader blockNorm structBook] = loadFileSBB(fileHeader, fileSBB);



if strcmp(type,'blocknorm')
    for k1=1:sbbHeader.numSignal
        figure, bar(blockNorm(k1,:))
        title(['Signal: ' num2str(k1)])
        xlabel('Block')
        ylabel('norm')
    end
end

min_amp = zeros(sbbHeader.numSignal,sbbHeader.numBlock);
max_amp = min_amp;
min_rho = min_amp;
max_rho = min_amp;
min_phase = min_amp;
max_phase = min_amp;

for k1=1:sbbHeader.numSignal
    for k2=1:sbbHeader.numBlock
        sb = structBook{k1,k2};
%         disp(['Signal: ' num2str(k1) '; Block: ' num2str(k2)])
%         disp(['min_amp: ' num2str(min(sb(:,2))) ]);
%         disp(['max_amp: ' num2str(max(sb(:,2))) ]);
%         disp(['min_rho: ' num2str(min(sb(:,3))) ]);
%         disp(['max_rho: ' num2str(max(sb(:,3))) ]);
%         disp(['min_phase: ' num2str(min(sb(:,5))) ]);
%         disp(['max_phase: ' num2str(max(sb(:,5))) ]);
        min_amp(k1,k2) = min(sb(:,2));
        max_amp(k1,k2) = max(sb(:,2));
        min_rho(k1,k2) = min(abs(sb(:,3)));
        max_rho(k1,k2) = max(abs(sb(:,3)));
        min_phase(k1,k2) = min(sb(:,5));
        max_phase(k1,k2) = max(sb(:,5));
    end
end

if strcmp(type,'parmrange')
    for k1=1:sbbHeader.numSignal
        figure; 
        
        subplot(3,1,1), bar([min_amp(k1,:)' max_amp(k1,:)'],'stacked')
        legend('min amp','max amp')
        title(['Signal: ' num2str(k1)])
        subplot(3,1,2), bar([min_rho(k1,:)' max_rho(k1,:)'],'stacked')
        legend('min rho','max rho')
        subplot(3,1,3), bar([min_phase(k1,:)' max_phase(k1,:)'],'stacked')
        legend('min phase','max phase')
        xlabel('Block')
        
    end
end

sb = structBook{signal,block};

if strcmp(type,'hist')
    

    if strcmp(parm,'coef')
        if isempty(dist)
            figure,hist(sb(:,2),nbin,dist);
        else
            figure,histfit(sb(:,2),nbin,dist);
        end
        title(['Signal: ' num2str(signal) '-- Block: ' num2str(block)])
    end

    if strcmp(parm,'decay')
        if isempty(dist)
            figure,hist(sb(:,3),nbin,dist);
        else
            figure,histfit(sb(:,3),nbin,dist);
        end
        title(['Signal: ' num2str(signal) '-- Block: ' num2str(block)])
    end

    if strcmp(parm,'freq')
        if isempty(dist)
            figure,hist(sb(:,4),nbin,dist);
        else
            figure,histfit(sb(:,4),nbin,dist);
        end
        title(['Signal: ' num2str(signal) '-- Block: ' num2str(block)])
    end

    if strcmp(parm,'phase')
        if isempty(dist)
            figure,hist(sb(:,5),nbin,dist);
        else
            figure,histfit(sb(:,5),nbin,dist);
        end
        title(['Signal: ' num2str(signal) '-- Block: ' num2str(block)])
    end

    if strcmp(parm,'timesupport')
        if isempty(dist)
            figure,hist(sb(:,6),nbin,dist);
        else
            figure,histfit(sb(:,6),nbin,dist);
        end
        title(['Signal: ' num2str(signal) '-- Block: ' num2str(block)])
    end
    
end


