function x = read_rdoptfile (fileName,numSignal) 


fid  = fopen(fileName,'r');


for k=1:numSignal
    string = fgets(fid); 
    string = fgets(fid);
    x{k} = sscanf(string,'%f');
end


fclose(fid);

