function formatFig(figHandler,figName,lang,sizeTick, sizeTitle, sizeLabel,figDim)
% This function was created by Leonardo Nunes (lonnes@lps.ufrj.br).
% This script configures the current figure.
% The folowing figure properties are connfigured:
% 
% text (ticks and labels) size
% text (ticks and labels) font
% line width
% figure size
% sizeTick, sizeTitle, sizeLabel are the sizes in pts for ticks, labels and title texts.

%
% if the variable 'figName' exists, a .fig and an .eps, files with the name
% contained in 'figName', are saved.

% Configuration:

font = 'Times New Roman';
lineWidth = 1;


%--------------------------------------------------------------------------
% Configuring figure

if(~isempty(get(0,'CurrentFigure')))
    
    set(figHandler,'Position',figDim);   
    fc = get(figHandler,'children'); % children of the current figure.

    % Cycling through children:
    
    for ii9183=1:length(fc)
        
        if(strcmp(get(fc(ii9183),'Type'),'axes'))

            % Configuring axes text:
            set(fc(ii9183),'FontSize',sizeTick);
            set(fc(ii9183),'FontName',font);
            
            % Configuring label text:
            ax = get(fc(ii9183),'xlabel');
            set(ax,'FontSize',sizeLabel);
            set(ax,'FontName',font);
            
            ay = get(fc(ii9183),'ylabel');
            set(ay,'FontSize',sizeLabel);
            set(ay,'FontName',font);       
            
            % Configuring title text:
            
            at = get(fc(ii9183),'title');
            set(at,'FontSize',sizeTitle);
            set(at,'FontName',font);
            
            ac = get(fc(ii9183),'children'); % axes children.
            
            for jj98719=1:length(ac)
               
                if(strcmp(get(ac(jj98719),'Type'),'line'))
                   
                    set(ac(jj98719),'LineWidth',lineWidth);
                    
                    if(strcmp(get(ac(jj98719),'marker'),'.'))
                        set(ac(jj98719),'markerSize',15);
                    end  
                    
                  
                end
                
                if(strcmp(get(ac(jj98719),'Type'),'text'))
                   
                    set(ac(jj98719),'FontSize',sizeLabel);
                    set(ac(jj98719),'FontName',font);
                    
                end
            end
            if(strcmpi(lang,'pt'))
                tick = get(fc(ii9183),'XTickLabel');
                tick = changeComma(tick);
                set(fc(ii9183),'XTickLabel',tick);
                tick = get(fc(ii9183),'YTickLabel');
                tick = changeComma(tick);
                set(fc(ii9183),'YTickLabel',tick);
                tick = get(fc(ii9183),'ZTickLabel');
                tick = changeComma(tick);
                set(fc(ii9183),'ZTickLabel',tick);
                set(fc(ii9183),'XTickMode','manual');
                set(fc(ii9183),'YTickMode','manual');
                set(fc(ii9183),'ZTickMode','manual');
            end
           
        end
    end
		set(figHandler,'PaperPositionMode','auto')
    set(figHandler,'Position',figDim);   
    aux = [figName '.eps'];
    saveas(figHandler,aux,'epsc2');
    aux = [figName '.fig'];

%  		pause(1)

    saveas(figHandler,aux,'fig');
        
end

function str = changeComma(str)
for ii = 1:size(str,1)

    str(ii,:) = regexprep(str(ii,:),'[.]', ',');

end
    
