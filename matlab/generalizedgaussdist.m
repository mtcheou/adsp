% Find the parameters of the generalized gaussian distribution of the
% vector parm
% mean, shape and scale


load decayparm.mat

% mean
mu = mean(parm);

% shape
p = [0.01:0.01:10];

Mp =  ( (gamma(2./p)).^2 ./ (gamma(1./p).*gamma(3./p)) ) ;


MpTarget = mean(abs(parm))^2 / mean(parm.^2);


ChosenShape =  p(find(abs(Mp-MpTarget)==min(abs(Mp-MpTarget))));

% scale

ChosenScale = sqrt( (var(parm)* gamma(1/ChosenShape)) / gamma(3/ChosenShape)  );


nbins =200;

x = -4:8/nbins:4;

xc = x(1:end-1)+4/nbins;

GGDx = (ChosenShape/(2*ChosenScale*gamma(1/ChosenShape)))*exp(- (abs(x-mu)/ ChosenScale).^(ChosenShape) );

GGD = (ChosenShape/(2*ChosenScale*gamma(1/ChosenShape)))*exp(- (abs(xc-mu)/ ChosenScale).^(ChosenShape) );


scaling = length(parm)/sum(GGD);


figure, hist(parm,nbins);
hold on;
plot(x,scaling*GGDx,'r--', 'LineWidth', 1.5);
hold off;


[N bin] = histc(parm,x);
figure, bar(xc,N(1:end-1),'grouped');
hold on;
plot(xc,scaling*GGD,'r--', 'LineWidth', 1.5);
hold off;


[h,p,st] = chi2gof(xc,'ctrs',xc,...
                  'frequency',N(1:end-1), ...
                  'expected',scaling*GGD, ...
                  'nparams',3);
