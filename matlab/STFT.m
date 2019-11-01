function [X PH] = STFT(x,ljanela,hop,nFFT)
%  X = STFT(x,ljanela,hop,nFrames,nFFT)
%  
%  x: signal to be decomposed
%  ljanela: Window length (in samples)
%  hop: Window advance
%  nFrames: number of frames
%  nFFT: size of FFT
%  
%  X: STFT
 
  w = hamming2(ljanela,hop);
  %  w = hamming(ljanela);
  T = length(x);
  nFrames = floor((T-ljanela)/hop+1);

  PH = zeros(nFFT,nFrames);
  X = zeros(nFFT,nFrames);
  X_aux = zeros(nFFT,nFrames);
  x = x(:);
  w = w(:);

  for m = 0:nFrames-1
    aux = [x(1+m*hop:m*hop+ljanela).*w; zeros(nFFT-ljanela,1)];
    aux = [aux(end/2+1:end); aux(1:end/2)];
    X_aux(:,m+1) = fft(aux,nFFT);
    X(:,m+1) = X_aux(:,m+1);
    PH(:,m+1) = phase(X_aux(:,m+1));
  end
 
 
 
end