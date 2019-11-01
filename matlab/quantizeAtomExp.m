function sbQ = quantizeAtomExp(step_amp,step_rho,step_phase,min_amp,min_rho,min_phase,sb)

sbQ(1) = linearquant(sb(1)-min_amp,step_amp) + min_amp;
sbQ(4) = linearquant(sb(4)-min_phase,step_phase) + min_phase;

if (sb(2)<0)
    sbQ(2) = - ( linearquant(abs(sb(2)) - min_rho,step_rho) + min_rho );
else
    sbQ(2) = linearquant(abs(sb(2))-min_rho,step_rho) + min_rho;
end

sbQ(3) = sb(3);
sbQ(5) = sb(5);
sbQ(6) = sb(6);
