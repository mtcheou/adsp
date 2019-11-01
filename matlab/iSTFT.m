function y = iSTFT(X, ljanela, window, hop, nFFT)
% y = DGistft(X, ljanela, window, hop, nFFT)
% Calculates the inverse STFT of X where a (ljanela)-point FFT has been calculated
% with window of ljanela points and hop number of points. 
%
% X = [B; B((end-1):-1:2,:)] where B is from specgram().
%
% David Gunawan 2004

y = zeros(size(X,2)*hop+ljanela-hop,1);
%  scale = mean(window.^2)*ljanela/hop;

H = ljanela/hop;

if(H == 2)
	gain_comp = zeros(ljanela/2,1);
	for k=1:ljanela/2
		for h = 0:H/2
			gain_comp(k) = gain_comp(k) + window(k+h*hop)^2;
		end
	end
else
	gain_comp = mean(window.^2)*ljanela/hop;
end
%  for k=1:ljanela/2
%      gain_comp(k) = 1/(window(k)^2+window(k+ljanela/2)^2);   
%  end


i = 1; % init index of start of frame
for f=1:size(X,2) % for every frame
   % ifft, window, overlap-add
%     y(i:i+ljanela-1) = y(i:i+ljanela-1) + real(ifft(X(:,f))).*window/scale; 
	 aux = real(ifft(X(:,f),nFFT));

	 aux = [aux(end/2+1:end); aux(1:end/2)];
    
   y(i:i+ljanela-1) = y(i:i+ljanela-1) + aux(1:ljanela).*window; 
	 y(i:i+ljanela/2-1) = y(i:i+ljanela/2-1) ./gain_comp;
   i = i + hop;
end
