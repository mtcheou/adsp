function playdecomp(matfile,FPS)

load(matfile)

%Define the area to be recorded 
hf=figure; 
rect = get(hf,'Position'); 
rect(1:2) = [0 0]; 

% Play the MATLAB movie 
clf 
N = 1;  
movie(hf,M,N,FPS,rect) 