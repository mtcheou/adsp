% audio_decomp.m

clear all
close all

N = 8;

dicParm = getDicParmDiscrete('std',N);

[numAtom numParm] = size(dicParm);
% figure;
% for ind = 1: numAtom
%     atomReal = genGaborReal(dicParm(ind,:),0 ,N);
%     plot(atomReal);
%     pause;
% end
chosenAtom = randperm(numAtom,1);
x = genGaborReal(dicParm(chosenAtom,:), 1.2 ,N);


% Fast MP Kolasa FFT
coef_vec = zeros(1,numAtom);
opt_phi_vec = zeros(1,numAtom);

deltaIndXi = N/2;
innerProd_xP = zeros(1,N/2);
innerProd_xQ = zeros(1,N/2);
innerProd_PQ = zeros(1,N/2);
normQ2 = zeros(1,N/2);
normP2 = zeros(1,N/2);
for ind = 1: deltaIndXi:numAtom
    s = dicParm(ind,1);
    u = dicParm(ind,1);
    
    g_window = genGaborReal([s u 0], 0,N);
    y1 = x.*g_window;
    Y1 = fft(y1);
    y2 = g_window.^2;
    Y2 = fft(y2);
    for k = 0 : N/2-1    
    %%%%%%%%%%%%% CONTINUAR AQUI
        xi = dicParm(ind,3);
        [opt_phi_vec(ind),coef_vec(ind)] = computeOptPhaseInnerProd(x,P,Q,xi);    
    end
end
ind_max = find(coef_vec==max(coef_vec));





% Fast MP FFT with time-shift




% Classic MP with optimum phase
% coef_vec = zeros(1,numAtom);
% opt_phi_vec = zeros(1,numAtom);
% for ind = 1: numAtom
%     atomComplex = genGaborComplex(dicParm(ind,:),N);
%     P = real(atomComplex);
%     Q = imag(atomComplex);    
%     innerProd_xP = sum(x.*P);
%     innerProd_xQ = sum(x.*Q);
%     innerProd_PQ = sum(P.*Q);
%     normQ2 = norm(Q)^2;
%     normP2 = norm(P)^2;    
%     xi = dicParm(ind,3);
%    [opt_phi_vec(ind),coef_vec(ind)] = computeOptPhaseInnerProd(innerProd_xP,innerProd_xQ,innerProd_PQ,normQ2,normP2,xi);    
% end
% ind_max = find(coef_vec==max(coef_vec));



for ind=1:length(ind_max)
    atomReal = genGaborReal(dicParm(ind_max(ind),:), opt_phi_vec(ind_max(ind)),N);
    figure,plot(0:N-1, [x' coef_vec(ind_max(ind))*atomReal']);
end







