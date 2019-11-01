function contAtomQuant = quantContAtomParm (contAtom)


nbits_amp =
max_amp =
min_amp =

nbits_rho =
max_rho =
min_rho =

nbits_phase =
max_phase =
min_phase =

nlevel_amp = (2.0^nbits_amp) - 1;
step_amp = ( max_amp - min_amp) / nlevel_amp ;

nlevel_rho = (2.0^nbits_rho) - 1;
step_rho =  (max_rho - min_rho) / nlevel_rho;

nlevel_phase = (2.0^nbits_phase) - 1;
step_phase = (max_phase - min_phase) / nlevel_phase;

for k1 =1:length(contAtom)
    sb = contAtom{k1};
    for k2 =1:length(sb)
        
    end
end