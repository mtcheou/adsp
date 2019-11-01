function plotrd(isignal,iblock, sigSize,fileName)

% plotrd - plots the rate-distortion operational curves
%
% Input:
%   isignal - signal number
%   iblock  - block number
%

marks = [ '.' 'o' 'x' '+' '*' 's' 'd' 'v' '^' '<' '>' 'p' 'h' '.' 'o' 'x' '+' '*' 's' 'd' 'v' '^' '<' '>' 'p'];

nmark = length(marks);

fname = [fileName '_rdpoint.out'];

[labels,x,y] = readColData(fname,11);

rdpoint = [x y];

% identifies rd points for each nb_amp
figure;
kk=1;
for k =4:2:12
    filter = find((rdpoint(:,1)==isignal)&(rdpoint(:,2)==iblock)&(rdpoint(:,3)==k)&(rdpoint(:,4)~=1)&(rdpoint(:,5)~=1));
    semilogy(rdpoint(filter,9)/sigSize,rdpoint(filter,10),marks(kk),'Color',[0 0 0]);
    kk=kk+1;
    hold on;
end
xlabel(' Taxa (bits por amostra)');
ylabel(' Erro quadrático (escala log)');
legend(cellstr(num2str([4:2:nmark]')));
hold off;

% nb_amp =8, vary rho
figure;
kk=1;
for k =4:2:22
    filter = find((rdpoint(:,1)==isignal)&(rdpoint(:,2)==iblock)&(rdpoint(:,3)==8)&(rdpoint(:,4)==k)&(rdpoint(:,5)~=1));
    semilogy(rdpoint(filter,9)/sigSize,rdpoint(filter,10),marks(kk),'Color',[0 0 0]);
    kk=kk+1;
    hold on;
end
xlabel(' Taxa (bits por amostra)');
ylabel(' Erro quadrático (escala log)');
legend(cellstr(num2str([4:2:nmark]')));
hold off;

% nb_amp =8, vary phi
figure;
kk=1;
for k =4:2:18
    filter = find((rdpoint(:,1)==isignal)&(rdpoint(:,2)==iblock)&(rdpoint(:,3)==8)&(rdpoint(:,4)~=1)&(rdpoint(:,5)==k));
    semilogy(rdpoint(filter,9)/sigSize,rdpoint(filter,10),marks(kk),'Color',[0 0 0]);
    kk=kk+1;
    hold on;
end
xlabel(' Taxa (bits por amostra)');
ylabel(' Erro quadrático (escala log)');
legend(cellstr(num2str([4:2:nmark]')));
hold off;


%%%%
filter = find((rdpoint(:,1)==isignal)&(rdpoint(:,2)==iblock)&(rdpoint(:,4)~=1)&(rdpoint(:,5)~=1));

figure,semilogy(rdpoint(filter,9)/sigSize,rdpoint(filter,10),'.');
hold on
xlabel(' Taxa (bits por amostra)');
ylabel(' Erro quadrático (escala log)');
%title(['Rate-distortion: ' 'signal: ' num2str(isignal) ' - block: ' num2str(iblock) ]);

fname = [fileName '_rdopcurve.out'];
[labels,x,y] = readColData(fname,12);

rdopcurve = [x y];

filter = find((rdopcurve(:,1)==isignal)&(rdopcurve(:,2)==iblock));

semilogy(rdopcurve(filter,9)/sigSize,rdopcurve(filter,10),'ro-');
hold off;

%%%%

% figure,semilogy(rdopcurve(filter,9)/sigSize,rdopcurve(filter,10)/(norma^2),'ro-');
% xlabel(' Taxa (bits por amostra)');
% ylabel(' Erro quadrático normalizado (escala log)');

