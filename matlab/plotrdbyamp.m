function plotrdbyamp(isignal,iblock,nbamp,iAmpRange,nbrho,nbphi,showby)

% plotrdbyamp - plots the rate-distortion operational curves of atoms
%               separated among amplitude ranges.
%
% Input:
%   isignal - signal number
%   iblock  - block number
%   nbamp   - number of bit for amplitude
%   iAmpRange - index of the amplitude range (grows with the amplitude level)
%   nbrho - chosen by encoder
%   nbphi - chosen by encoder
%   showby - 'vary_phase', 'vary_rho'
%

if strcmp(showby,'vary_rho')
    chosenparm = 5;
elseif strcmp(showby,'vary_phase')
    chosenparm = 6;
end


[labels,x,y] = readColData('rdbyamp_rdpoint.out',10,0);

rdpoint = [x y];


x=[];
y=[];
for k =1:16

    filter = find(  (rdpoint(:,1)==isignal)& ...
                    (rdpoint(:,2)==iblock) & ...
                    (rdpoint(:,3)==nbamp) &...
                    (rdpoint(:,4)==iAmpRange) & ...
                    (rdpoint(:,chosenparm)==k));

    x = [x rdpoint(filter,9) ];
    y = [y rdpoint(filter,10)];
end           
            
%rdpointfilt = rdpoint(filter,:);


            

%figure,plot(rdpoint(filter,9),rdpoint(filter,10),'*');
figure,plot(x,y,'*');
hold on
xlabel(' Rate (total bit number)');
ylabel(' Distortion (MSE)');
title(['Rate-distortion: ' 'signal: ' num2str(isignal) ... 
       ' - block: ' num2str(iblock) ...
       ' - nbamp: ' num2str(nbamp) ...
       ' - iAmpRange: ' num2str(iAmpRange) ]);
   
%legend('8','10','12','14','16')

[labels,x,y] = readColData('rdbyamp_rdopcurve.out',11,0);

rdopcurve = [x y];

filter = find(  (rdopcurve(:,1)==isignal)& ...
                (rdopcurve(:,2)==iblock) & ...
                (rdopcurve(:,3)==nbamp) &...
                (rdopcurve(:,4)==iAmpRange));
            
rdopcurvefilt = rdopcurve(filter,:);

plot(rdopcurve(filter,9),rdopcurve(filter,10),'ro-');

if (nargin>4)
    filter = find(  (rdpoint(:,1)==isignal)& ...
                    (rdpoint(:,2)==iblock) & ...
                    (rdpoint(:,3)==nbamp) &...
                    (rdpoint(:,4)==iAmpRange) & ...
                    (rdpoint(:,5)==nbrho) & ...
                    (rdpoint(:,6)==nbphi) );
    plot(rdpoint(filter,9),rdpoint(filter,10),'k<','MarkerSize',12);
end



hold off;