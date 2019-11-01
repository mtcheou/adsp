function y = linearquant(x,step)

if (step==0.0)
    y = x;    
else
    y = floor( (x + step/2) / step) * step;
end
