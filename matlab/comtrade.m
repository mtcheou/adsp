function [v,f] = comtrade( arquivo, Fs, dados)
% function comtrade( arquivo, Fs, dados)
%    gera arquivo COMTRADE ASCII a partir de uma matriz de dados
%      onde cada linha e' um vetor distinto
% function v = comtrade(arquivo)
%    le arquivo comtrade ASCII e coloca na matrix v. O valor da freq. de amostragem 
%      vai em f

% Autor : Marco Antonio M. Rodrigues
% novembro/1995

f = 0;
v = 0;

if nargout<1	% está escrevendo no arquivo

	% trata os dados
	[L C] = size(dados);
	if L>C
		dados = dados';
	end
	a = max([L C]);
	L = min([L C]);
	C = a;
	ganho = max(dados');
	minimo = min(dados');
	ganho = max(abs(ganho),abs(minimo));
	[a,t]=min(ganho);
	while a==0,
		ganho(t)=1;
		[a,t]=min(ganho);
	end
	ganho = ganho/32767;
	dados = diag(1./ganho) * dados;
	% campos do comtrade
	a=(1:C);
	t=(0:C-1)*1E6/Fs;
	string = '%010.0f,%010.0f';
	dad = [a;t];
	for i=1:L
	   string = [string ',%06.0f'];
	   dad = [dad;dados(i,:)];
	end
	string = [string '\r\n'];
	% escreve no arquivo de dados
	a=findstr(arquivo,'.');
	if ~isempty(a)
	   arquivo = [arquivo(1:(a-1)) '.dat'];
	else
	   arquivo = [arquivo '.dat'];
	end
	fid = fopen(arquivo,'w');
	fprintf( fid, string, dad);
	fprintf( fid,'\r\n');
	fclose(fid);

	% escreve no arquivo de configuracao
	a=findstr(arquivo,'.');
	arquivo = [arquivo(1:(a-1)) '.cfg'];
	fid = fopen(arquivo,'w');
	fprintf( fid,'%s,MATLAB\r\n', arquivo(1:(a-1))); % substacao,RDP
	t = int2str(L);
	fprintf( fid,'%d,%dA,0D\r\n', L, L); % num canais
	for i=1:L,   % canais
		fprintf(fid,'%d,nome%d,,,V,%G,0,0,-32767,32767\r\n',i,i,ganho(i));
	end
	fprintf( fid,'60\r\n'); % freq da rede
	fprintf( fid,'1\r\n'); % numero de taxas
	fprintf( fid,'%G,%d\r\n', Fs, C);  % taxas
	fprintf( fid,'11/11/95,12:00:00.000000\r\n');
	fprintf( fid,'11/11/95,12:00:00.000000\r\n');
	fprintf( fid,'ASCII\r\n\r\n');
	fclose(fid);

elseif nargout>=1	% está lendo do arquivo

	% le arquivo de configuração
	a=findstr(arquivo,'.');
	if ~isempty(a)
	   arquivo = [arquivo(1:(a-1)) '.cfg'];
	else
	   arquivo = [arquivo '.cfg'];
	end
	fid = fopen(arquivo,'r');
	% Le primeira linha
	texto = fgets(fid);
	%disp(texto)
	% Le número de canais
	dummy = fgets(fid);	a=[findstr(dummy,',') length(dummy)-2];
	TotCan = str2num(dummy(1:a(1)-1));
	AnaCan = str2num(dummy(a(1)+1:a(2)-2));
	DigCan = str2num(dummy(a(2)+1:a(3)-1));
	% Lê características dos canais
	for i=1:TotCan
		dummy=fgets(fid);	a=[findstr(dummy,',') length(dummy)-2];
		canal(i) = str2num(dummy(1:a(1)-1));
		if i<=AnaCan
			cang(i) = str2num(dummy(a(5)+1:a(6)-1));
			clin(i) = str2num(dummy(a(6)+1:a(7)-1));
		end		
	end
	FreqRede = str2num(fgets(fid));
	NFreq = str2num(fgets(fid));
	if NFreq~=1
	    error('COMTRADE: este programa não permite variar a freq. de amostragem')
	  end
	dummy=fgets(fid);
	a=findstr(dummy,',');
	if ~isempty(a)
	   Fs = str2num(dummy(1:(a-1)));
	   NAmo = str2num(dummy(a+1:end));
	else
	   error('COMTRADE: erro de formato no campo de freq. de amostragem')
	end
	for i=1:2, dummy=fgets(fid);, end
	dummy=fgetl(fid);
	if strcmp(dummy,'ASCII')
	  tipo = 'A';
	elseif strcmp(dummy,'BINARY')
	  tipo = 'B';
	else
	  error('COMTRADE: erro de formato no campo de tipo de arquivo de dados')
	end
	fclose(fid);
	a=findstr(arquivo,'.');
	if ~isempty(a)
	   arquivo = [arquivo(1:(a-1)) '.dat'];
	else
	   arquivo = [arquivo '.dat'];
	end
	fid = fopen(arquivo,'r');
	v =[];
	if (tipo=='A')
		for i=1:NAmo
			dummy = fgets(fid);
			a=[findstr(dummy,',') length(dummy)-1];
			for j=1:TotCan
		    	v(i,j) = str2num(dummy(a(1+j)+1:a(2+j)-1));
				if j<=AnaCan
					v(i,j) = v(i,j)  * cang(j) + clin(j);
				end
			end
		end
	end
	
	if (tipo=='B')
		s = fread(fid,'uint16')
		d = downsample(s,5,4);
		figure,plot(d);
		for j=1:AnaCan
			v =[v;  d(j+4,:)*cang(j)+clin(j)]; 
		end
	end
	%cang;
	%clin;
	fclose(fid);
	f = Fs;
end