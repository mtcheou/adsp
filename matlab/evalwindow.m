

ljanela = 200;
hop = 100;

%w = hamming2(ljanela,hop);
w = hamming(ljanela);
w = hann(ljanela);

T = 1000;
nFrames = floor((T-ljanela)/hop+1);


y =zeros(1,T);

figure;
for m = 0:nFrames-1
    plot(m*hop +1 : m*hop +ljanela,w,'o');
    y(m*hop +1 : m*hop +ljanela) = y(m*hop +1 : m*hop +ljanela) + w';
    hold on;
end

hold off;
 
figure, plot(y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 30;          % window length
R = (M)/3;     % hop size
N = 3*M;         % overlap-add span
w = blackman(M);  % window
z = zeros(N,1);  plot(z,'-k');  hold on;  s = z;
for so=0:R:N-M
  ndx = so+1:so+M;        % current window location
  s(ndx) = s(ndx) + w;    % window overlap-add
  wzp = z; wzp(ndx) = w;  % for plot only 
  plot(wzp,'--ok');       % plot just this window
end
plot(s,'ok');  hold off;  % plot window overlap-add
