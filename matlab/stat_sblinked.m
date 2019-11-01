load sblinked_load.out

nb_coef = 8;
nb_rho = 12;
nb_phi = 12;
figure;
%%%
coef = sblinked_load(:,1);
coef_max = max(coef);
delta_coef = coef_max/((2^nb_coef)-1);
t = nb_coef:-1:0;
edges = [0 coef_max*(2.^(-t))];
n = histc(coef,edges);
% 
bar(edges,n,'histc');
title('Coef');

%histfit(coef,100,'gamma')


nel = length(edges)-1;
for k = 2:nel-1
    ind_range{k-1} = find(coef>=edges(k-1)&coef<edges(k));
end
ind_range{nel-1} = find(coef>=edges(nel-1)&coef<=edges(nel));

coef_edges = edges;
pause;

%%%
rho = abs(sblinked_load(:,2));
rho_max = max(rho);
delta_rho = rho_max/((2^nb_rho)-1);
t = 0:(2^nb_rho)-1;
edges = t*delta_rho;
n = histc(rho,edges);
 bar(edges,n,'histc');
title('Abs(Rho)')

rho_edges = edges;
pause;


%%%
freq_tab = tabulate(sblinked_load(:,3));
 bar((44100/(2*pi))*freq_tab(:,1),freq_tab(:,3))
xlabel('Freq (Hz)')
ylabel('%')
title('Freq')
pause;

% figure, bar((44100/(2*pi))*freq_tab(:,1),freq_tab(:,2))
% xlabel('Freq (Hz)')
% ylabel('Count')
% pause;

%%%
phi = abs(sblinked_load(:,4));
phi_max = max(phi);
delta_phi = phi_max/((2^nb_phi)-1);
t = 0:(2^nb_phi)-1;
edges = t*delta_phi;
n = histc(phi,edges);
 bar(edges,n,'histc');
title('Phase')
phi_edges=edges;
pause;

%%%
dt_tab = tabulate(sblinked_load(:,6)-sblinked_load(:,5));
 bar(dt_tab(:,1),dt_tab(:,3))
xlabel('Delta t (samples)')
ylabel('%')
title('Delta t')
pause;

% figure, bar(dt_tab(:,1),dt_tab(:,2))
% xlabel('Delta t (samples)')
% ylabel('Count')
% title('Delta t')
% pause;

%%%
figure, scatter(sblinked_load(:,6)-sblinked_load(:,5),abs(sblinked_load(:,2)))
xlabel('Delta t')
ylabel('Abs(Rho)')
pause;

%%%
nrange = nel - 1;
for k =1:nrange
    n = histc(rho(ind_range{k}),rho_edges);
     bar(rho_edges,n,'histc');
    title(['Rho - coef range: ' num2str(k) ': ' num2str(coef_edges(k)) ' - ' num2str(coef_edges(k+1)) ] )
    pause;
end

for k =1:nrange
    freq_tab = tabulate(sblinked_load(ind_range{k},3));
     bar((44100/(2*pi))*freq_tab(:,1),freq_tab(:,3))
    xlabel('Freq (Hz)')
    ylabel('%')
    title(['Freq - coef range: ' num2str(k) ': ' num2str(coef_edges(k)) ' - ' num2str(coef_edges(k+1)) ] )
    pause;
end

for k =1:nrange
    n = histc(phi(ind_range{k}),phi_edges);
     bar(phi_edges,n,'histc');
    title(['Phi - coef range: ' num2str(k) ': ' num2str(coef_edges(k)) ' - ' num2str(coef_edges(k+1)) ] )
    pause;
end

for k =1:nrange
    dt_tab = tabulate(sblinked_load(ind_range{k},6)-sblinked_load(ind_range{k},5));
     bar(dt_tab(:,1),dt_tab(:,3))
    xlabel('Delta t (samples)')
    ylabel('%')
    title(['Delta t - coef range: ' num2str(k) ': ' num2str(coef_edges(k)) ' - ' num2str(coef_edges(k+1)) ] )
    pause;
end


%%%%

dt = sblinked_load(:,6)-sblinked_load(:,5);
hist3([rho dt],[100 100])

