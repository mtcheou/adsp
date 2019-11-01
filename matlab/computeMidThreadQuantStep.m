function step =  computeMidThreadQuantStep(nb,min,max)

nlevel = ((2.0)^nb)  - 1.0;
if (nlevel==0.0)
    step = 1e+30;
else
    step = abs( (max - min) / nlevel ); 
end
