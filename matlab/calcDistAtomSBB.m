function matcell = calcDistAtomSBB(fileHeader, fileSBB, nbamp,nbrho,nbphi )


[sbbHeader blockNorm structBook] = loadFileSBB(fileHeader, fileSBB);

sbbHeader.blockSize

for k1=1:sbbHeader.numSignal
    for k2=1:sbbHeader.numBlock
        sb = structBook{k1,k2};
        sbQ = sb;
        min_amp = min(sb(:,2));
        max_amp = max(sb(:,2));
        min_rho = min(abs(sb(:,3)));
        max_rho = max(abs(sb(:,3)));
        min_phi = min(sb(:,5));
        max_phi = max(sb(:,5));
        [r c] = size(sb)
        
        step_amp = computeMidThreadQuantStep(nbamp,min_amp,max_amp);
        step_rho = computeMidThreadQuantStep(nbrho,min_rho,max_rho);
        step_phi = computeMidThreadQuantStep(nbphi,min_phi,max_phi);
        
        mat = zeros(r,(1+6+6+6));
        for k3 =1:r
            sbQ(k3,2:end) = quantizeAtomExp(step_amp,step_rho,step_phi, ...
                                            min_amp,min_rho,min_phi, ...
                                            sb(k3,2:end));
            [df kappa2 kappaq2] = calcAtomDistFunction(sb(k3,2:end), ...
                                                        sbQ(k3,2:end), ...
                                                        sbbHeader.blockSize);
            [x x_norm] = genexp(sb(k3,3:end),sbbHeader.blockSize);
            
            [x_q x_q_norm] = genexp(sbQ(k3,3:end),sbbHeader.blockSize);
            
            %distAtom = sum((x_norm*x - x_q_norm*x_q).^2 )/(x_norm*x_norm);
            distAtom = norm(x-x_q)^2;
            
            mat(k3,:) = [k3 sb(k3,2:end) sbQ(k3,2:end) distAtom df (x_norm^2) kappa2 (x_q_norm^2) kappaq2];
        end
        structBookQ{k1,k2} = sbQ;
        matcell{k1,k2} = mat;
    end
end


