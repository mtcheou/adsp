function [bloco_onset,onsets,hop,num_blocos] = onsetsf(x)
% [bloco_onset,onsets,HOP] = bloco_onsetF(x);
% Escrita por Alexandre Leizor em dezembro de 2007.
% Função para detecção de bloco_onset de sinais musicais.
% De acordo com a técnica spectral flux do artigo Onset Detection Revisited do Simon Dixon DAFx-2006.
% x = sinal de entrada (vetor coluna)

N = 2048; % Como no artigo (pode virar um dos argumentos da função)
hop = 441; % idem
comp = length(x);
excesso = mod(comp-N,hop); % Comprimento do trecho de sinal ao fim de x que não completa o último bloco com N amostras.
sinal = x(1:comp-excesso); % Truncamento do sinal para evitar que a análise exceda o comprimento do sinal. 
num_blocos = ((length(sinal)-N)/hop)+1;
tamanho_da_DFT = N;

% No artigo os cálculos são realizados usando todo o espectro, mas como o sinal é real,...
% ... só metade precisa ser usada já que o espectro é simétrico.
DFTs = zeros(tamanho_da_DFT/2-1,num_blocos); % Cada coluna receberá uma transformada de N bins.

for bloco = 1:(num_blocos),
    fft_temp = fft(sinal((bloco-1)*hop+1:(bloco-1)*hop+N).*hamming(N),tamanho_da_DFT);
    DFTs(:,bloco) = abs(fft_temp(2:tamanho_da_DFT/2));
end
clear fft_temp

argumento = diff(DFTs,1,2); % Argumento da função retificadora.
H = (argumento + abs(argumento))/2;
SF = sum(H,1);
SF = SF - mean(SF);
SF = SF/std(SF);

% Suavização com média móvel
m = 3;
w = 3;
SF_suavizado = zeros(length(SF));
SF = [SF(m*w+1:-1:2) SF SF(length(SF)-1:-1:length(SF)-w)]; % Acrescenta m*w amostras espelhadas em torno da 1a amostra no início e w amostras espelhadas em torno da amostra final do sinal para permitir a suavização em toda a extensão de SF.
for c = 1:length(SF_suavizado);
    c_sf = c+m*w;
    SF_suavizado(c) = mean(SF(c_sf-m*w:c_sf+w));
end
SF = SF(m*w+1:length(SF)-w);
delta = 1;
limiares = delta + SF_suavizado;

onsets = zeros(1,length(x));
% O for abaixo serve para implementar a primeira condição do item 2.6 do artigo do Dixon
for n = 1:length(SF),
    if (n < w+1),
       if ((SF(n)>limiares(n)) && isempty(find(SF(1:n-1)>SF(n))) && isempty(find(SF(n+1:n+w)>SF(n)))),
          onsets(floor(n*hop)) = 1;
       end 
    elseif (n > length(SF)-w),
       if ((SF(n)>limiares(n)) && isempty(find(SF(n-w:n-1)>SF(n))) && isempty(find(SF(n+1:length(SF))>SF(n)))),
          onsets(floor(n*hop)) = 1;
       end
    else
       if ((SF(n)>limiares(n)) && isempty(find(SF(n-w:n-1)>SF(n))) && isempty(find(SF(n+1:n+w)>SF(n)))),
          onsets(floor(n*hop)) = 1;
       end
    end
end

bloco_onset = zeros(1,num_blocos);
for bloco = 1:num_blocos,
    if onsets(floor(bloco*hop)) == 1;
       bloco_onset(bloco) = 1;
    end
end

% A condição 3, opcional, não foi implementada.

figure
subplot(2,1,1)
plot(sinal)
axis('tight')
hold on
stem(.5*onsets(1:length(sinal)),'--k')
hold off
subplot(2,1,2)
plot(SF)
axis('tight')
hold on
plot(limiares,'g');
hold off
